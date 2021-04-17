import numpy as np
import math
from datetime import datetime
import collections
from operator import itemgetter
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd

from Network import Network


def main(machine_input_file, travel_time_file):

    # Simulation parameters
    max_time = 100000
    warm_up_time = 0
    nr_runs = 1
    policy = "GREEDY"  # Set to "REACTIVE" or "GREEDY"
    dependency = True  # Set to False for Question 1, True for Question 2

    # Read the input files
    nr_machines, list_thresholds, list_lambdas, list_alphas, list_betas, list_preventative_times, \
        list_corrective_times, list_preventative_costs, list_corrective_costs, list_downtime_costs, \
        list_travel_times = read_input(machine_input_file, travel_time_file)

    # Initializing the average variables storing data for each run
    average_operational_percentages = []
    average_response_times = []
    average_maintenance_costs = []
    average_total_costs = []

    for i in range(nr_runs):

        print("RUN %s" % (i + 1))

        # Create the network
        network = Network(nr_machines, policy, list_thresholds, list_lambdas, list_alphas, list_betas,
                          list_preventative_times, list_corrective_times, list_preventative_costs,
                          list_corrective_costs, list_downtime_costs, list_travel_times, dependency, warm_up_time)

        # Simulate network and get output variables
        non_operational_time, response_times, maintenance_costs, downtime_costs, non_operational_times_stamped, \
            maintenance_costs_stamped = simulate_network(network, max_time)

        # Calculating average operational percentages for the run
        average_operational_percentages.append([1 - (i / (max_time - warm_up_time)) for i in non_operational_time])

        # Calculating average response times for the run
        average_response_times.append([sum(i) / len(i) for i in response_times])

        # Calculating average maintenance costs for the run
        average_maintenance_costs.append([i / (max_time - warm_up_time) for i in np.add(maintenance_costs, downtime_costs)])

        # Calculating total costs for the run
        average_total_costs.append(np.sum(np.add(maintenance_costs, downtime_costs)) / (max_time - warm_up_time))

        # Histogram of Response Times for each machine
        fig, axs = plt.subplots(network.nr_machines, figsize=(15, 15))
        fig.suptitle('Response Times for %s Policy, Dependency = %s ' % (policy, dependency))
        for machine in range(network.nr_machines):

            axs[machine].hist(response_times[machine], bins=20, ec='black')
            axs[machine].set_title("Machine %s" % (machine))

        fig.tight_layout()
        fig.subplots_adjust(top=0.93)
        plt.show()
        fig.savefig('response_times_%s_%s.png' % (policy, dependency), bbox_inches='tight', dpi=300)

        # Plotting Average Mainteance Costs using time stamped output
        fig_maintenance, axs = plt.subplots(network.nr_machines, figsize=(15, 15))
        fig_maintenance.suptitle('Average Mainteance Costs vs Time for %s Policy, Dependency = %s ' % (policy, dependency))

        for machine in range(network.nr_machines):

            # Specifying x and y values to plot
            x = [i[1] for i in maintenance_costs_stamped[machine]]
            y = [i[0] / i[1] for i in maintenance_costs_stamped[machine]]

            # Plotting and formatting the plot
            axs[machine].plot(x, y)
            axs[machine].set_title("Machine %s" % (machine))
            axs[machine].xaxis.set_major_locator(ticker.MultipleLocator(5000))

        fig_maintenance.subplots_adjust(top=0.93)
        plt.show()
        fig_maintenance.savefig('maintenance_costs_%s_%s.png' % (policy, dependency), bbox_inches='tight', dpi=300)

        # Plotting Average Non Operational Percentages using time stamped output
        fig_non_operational, axs = plt.subplots(network.nr_machines, figsize=(15, 15))
        fig_non_operational.suptitle('Average Non Operational Percentage vs Time for %s Policy, Dependency = %s ' % (policy, dependency))

        for machine in range(network.nr_machines):

            # Specifying x and y values to plot
            x = [i[1] for i in non_operational_times_stamped[machine]]
            y = [i[0] / i[1] for i in non_operational_times_stamped[machine]]

            # Plotting and formatting the plot
            axs[machine].plot(x, y)
            axs[machine].set_title("Machine %s" % (machine))
            axs[machine].xaxis.set_major_locator(ticker.MultipleLocator(5000))

        fig_non_operational.subplots_adjust(top=0.93)
        plt.show()
        fig_non_operational.savefig('non_operational_%s_%s.png' % (policy, dependency), bbox_inches='tight', dpi=300)

        print("Average Operational Times:", [1 - (i / (max_time - warm_up_time)) for i in non_operational_time])
        print("Average Response Times:", [sum(i) / len(i) for i in response_times])
        print("Average Maintenance Costs:", [i / (max_time - warm_up_time) for i in maintenance_costs])
        print("Average Downtime Costs:", [i / (max_time - warm_up_time) for i in downtime_costs])
        print()

        # Deleting the network before the next run
        del network

    # Initializing output table that is going to be used for exporting the results to excel file
    output_table_dictionary = {}
    standard_deviations = []
    for machine in range(nr_machines):

        # Calculating the averages for output variables per machine using all runs
        average_operational_percentages_for_station = [item[machine] for item in average_operational_percentages]
        average_response_times_for_station = [item[machine] for item in average_response_times]
        average_maintenance_costs_for_station = [item[machine] for item in average_maintenance_costs]

        # Calculating the standard deviation for output variables using all runs
        std_operational_percentages_for_station = round(np.std(average_operational_percentages_for_station), 4)
        std_response_times_for_station = round(np.std(average_response_times_for_station), 4)
        std_maintenance_costs_for_station = round(np.std(average_maintenance_costs_for_station), 4)
        std_total_costs = round(np.std(average_total_costs), 4)

        # Calculating the variance for output variables using all runs
        variance_operational_percentages_for_station = round(std_operational_percentages_for_station ** 2, 4)
        variance_response_times_for_station = round(std_response_times_for_station ** 2, 4)
        variance_maintenance_costs_for_station = round(std_maintenance_costs_for_station ** 2, 4)
        variance_total_costs = round(std_total_costs ** 2, 4)

        # Adding standard deviations to use it for estimating standard deviation on the required number of runs equation

        standard_deviations.extend([std_operational_percentages_for_station, std_response_times_for_station,
                                    std_maintenance_costs_for_station])

        # Calculating confidence intervals
        operational_percentages_lower, operational_percentages_upper, operational_percentages_mean = \
            calc_conf_interval_list(average_operational_percentages_for_station, nr_runs)
        response_times_lower, response_times_upper, response_times_mean = \
            calc_conf_interval_list(average_response_times_for_station, nr_runs)
        maintenance_costs_lower, maintenance_costs_upper, maintenance_costs_mean = \
            calc_conf_interval_list(average_maintenance_costs_for_station, nr_runs)
        total_costs_lower, total_costs_upper, total_costs_mean = calc_conf_interval_list(average_total_costs, nr_runs)

        print(f"Average Operational Percentage for machine {machine} is {operational_percentages_mean} on average,standard_deviation = {std_operational_percentages_for_station}, lower = {operational_percentages_lower}, upper = {operational_percentages_upper}")
        print(f"Average Response Time for machine {machine} is {response_times_mean} on average, standard_deviation = {std_response_times_for_station}, lower = {response_times_lower}, upper = {response_times_upper}")
        print(f"Average Maintenance Cost for machine {machine} is {maintenance_costs_mean} on average, standard_deviation = {std_maintenance_costs_for_station}, lower = {maintenance_costs_lower}, upper = {maintenance_costs_upper}")
        print()

        output_table_dictionary[machine] = [operational_percentages_lower, operational_percentages_mean,
                                            operational_percentages_upper, variance_operational_percentages_for_station,
                                            response_times_lower, response_times_mean,
                                            response_times_upper, variance_response_times_for_station,
                                            maintenance_costs_lower, maintenance_costs_mean,
                                            maintenance_costs_upper, variance_maintenance_costs_for_station,
                                            total_costs_lower, total_costs_mean, total_costs_upper, variance_total_costs
                                            ]

    # Convert the results dictionary to dataframe and export it to excel file
    df = pd.DataFrame.from_dict(output_table_dictionary, orient='index',
            columns=['op-lower', 'op-mean', 'op-upper', 'op-variance',
                    'r-lower', 'r-mean', 'r-upper', 'r-variance',
                    'ac-lower', 'ac-mean', 'ac-upper', 'ac-variance',
                    'tc-lower', 'tc-mean', 'tc-upper', 'tc-variance'])

    df.to_excel("output %s_%s.xlsx" % (policy, dependency))


# Method to simulate the network with given parameters
def simulate_network(network, max_time):

    while network.current_time < max_time:

        # Handling the first decision point(t = 0)
        if network.current_time == 0:

            if network.policy == "REACTIVE":

                selected_machine, failed_machines_count = network.find_max_downtime()

                if failed_machines_count == 0:

                    # Since no machine is failed at the start FSE decides to stay IDLE
                    network.fse.state = "IDLE"

                    # Updating the last idle time statistic
                    network.last_idle_time = network.current_time

                else:

                    print("Error in initializing the network")

            elif network.policy == "GREEDY":

                selected_machine, functional_or_not = network.find_biggest_degradation()

                if selected_machine == -1:

                    # Since no machine is non-healthy at the start fse decides to stay IDLE
                    network.fse.state = "IDLE"

                else:

                    print("Error in initializing the network")

        # Get earliest event from the front of the FES
        current_event = network.fes.next()

        # Update time based on current_event
        network.current_time = current_event.time

        if network.policy == "REACTIVE":

            # Checking if the current event is DEGRADATION
            if current_event.typ == "DEGRADATION":

                # Handle the degradation
                network.handle_degradation(current_event)

                # Handling the failure case for the machine
                if network.machines[current_event.machine].degradation >= network.machines[current_event.machine].failure_threshold:

                    # If degradation exceeds failure threshold set it to the failure threshold
                    network.machines[current_event.machine].degradation = network.machines[current_event.machine].failure_threshold

                    # Setting the state of the machine to failed
                    network.machines[current_event.machine].state = "FAILED"

                    # If there is dependency update all the arrival intensities
                    if network.dependency is True:

                        network.handle_dependency(network.lambdas)

                    # Updating the last_non_operational_times array
                    network.last_non_operational_times[current_event.machine] = network.current_time

                    # Checking if the field service engineer is busy
                    if network.fse.state == "IDLE":

                        # Handle fse action for the current event
                        network.handle_fse_action_reactive(current_event)

                # If failure doesn't occur schedule the next degradation event for the machine
                else:

                    network.schedule_degradation(current_event.machine)

            # Checking if the current event is ARRIVAL
            elif current_event.typ == "ARRIVAL":

                # Update the total_travel_time statistic if the simulation time is above warm-up time
                if network.current_time > network.warm_up_time:

                    # Updating total_travel_time statistic
                    network.total_travel_time += network.current_time - network.last_travel_time

                # If arrival event occurs update the location of the field service engineer
                network.fse.location = current_event.machine

                # Handle fse action for the current event
                network.handle_fse_action_reactive(current_event)

            # Checking if the current event is END_MAINTENANCE
            elif current_event.typ == "END_MAINTENANCE":

                # Update the performance measures if the simulation time is above warm-up time
                if network.current_time > network.warm_up_time:

                    # Updating the non_operational_times statistic
                    if network.last_non_operational_times[current_event.machine] > network.warm_up_time:

                        network.non_operational_times[current_event.machine] += network.current_time - network.last_non_operational_times[current_event.machine]

                    else:

                        network.non_operational_times[current_event.machine] += network.current_time - network.warm_up_time

                    # Update time stamped non_operational_time statistic
                    network.non_operational_times_stamped[current_event.machine].append((network.non_operational_times[current_event.machine], network.current_time))

                    # Updating the number_of_repairs statistic
                    network.number_of_repairs[current_event.machine] += 1

                    # Updating the total downtime cost
                    downtime_cost_machine = network.downtime_costs[current_event.machine]
                    network.total_downtime_costs[current_event.machine] += downtime_cost_machine * math.floor(network.current_time - network.last_non_operational_times[current_event.machine])

                    # Update time stamped maintenance cost statistic
                    network.maintenance_costs_stamped[current_event.machine].append((network.maintenance_costs[current_event.machine] + network.total_downtime_costs[current_event.machine], network.current_time))

                    # Updating total_repair_time statistic
                    network.total_repair_time += network.current_time - network.last_repair_time

                # Setting the current machine's state to "FUNCTIONAL" and resetting its degradation level
                network.machines[current_event.machine].state = "FUNCTIONAL"
                network.machines[current_event.machine].degradation = 0

                # If there is dependency update all the arrival intensities
                if network.dependency is True:

                    network.handle_dependency(network.lambdas)

                # Schedule the next degradation event for the machine
                network.schedule_degradation(current_event.machine)

                # Handle fse action for the current event
                network.handle_fse_action_reactive(current_event)

        if network.policy == "GREEDY":

            if current_event.typ == "DEGRADATION" and network.machines[current_event.machine].degradation == network.machines[current_event.machine].failure_threshold:

                network.removed_event_count += 1

            # If the current event is degradation and the machine is not being repaired or failed at the moment
            if current_event.typ == "DEGRADATION" and network.machines[current_event.machine].degradation != network.machines[current_event.machine].failure_threshold:

                if network.fse.location == current_event.machine and network.fse.state == "REPAIR":

                    print("BIG ERROR")

                # Handle the degradation
                network.handle_degradation(current_event)

                # Handling the failure case for the machine
                if network.machines[current_event.machine].degradation >= network.machines[current_event.machine].failure_threshold:

                    # If degradation exceeds failure threshold it is now equal to the threshold
                    network.machines[current_event.machine].degradation = network.machines[current_event.machine].failure_threshold

                    # Setting the state of the machine to failed
                    network.machines[current_event.machine].state = "FAILED"

                    # If there is dependency update all the arrival intensities
                    if network.dependency is True:

                        network.handle_dependency(network.lambdas)

                    # Updating the last_non_operational_times array
                    network.last_non_operational_times[current_event.machine] = network.current_time

                # If failure doesn't occur schedule the next degradation event for the machine
                else:

                    network.schedule_degradation(current_event.machine)

                # Checking if the field service engineer is busy
                if network.fse.state == "IDLE":

                    # If the fse is idle handle the fse action
                    network.handle_fse_action_greedy(current_event)

            elif current_event.typ == "ARRIVAL":

                # Update the total_travel_time statistic if the simulation time is above warm-up time
                if network.current_time > network.warm_up_time:

                    # Updating total_travel_time statistic
                    network.total_travel_time += network.current_time - network.last_travel_time

                # If arrival event occurs update the location of the field service engineer
                network.fse.location = current_event.machine

                # Handle the fse action
                network.handle_fse_action_greedy(current_event)

            elif current_event.typ == "END_MAINTENANCE":

                # Update the performance measures if the simulation time is above warm-up time
                if network.current_time > network.warm_up_time:

                    # Updating the non_operational_times statistic
                    if network.last_non_operational_times[current_event.machine] > network.warm_up_time:

                        network.non_operational_times[current_event.machine] += network.current_time - network.last_non_operational_times[current_event.machine]

                    else:

                        network.non_operational_times[current_event.machine] += network.current_time - network.warm_up_time

                    # Update time stamped non_operational_time statistic
                    network.non_operational_times_stamped[current_event.machine].append((network.non_operational_times[current_event.machine], network.current_time))

                    # Updating the number_of_repairs statistic
                    network.number_of_repairs[current_event.machine] += 1

                    # Updating the total downtime cost
                    downtime_cost_machine = network.downtime_costs[current_event.machine]
                    network.total_downtime_costs[current_event.machine] += downtime_cost_machine * math.floor(network.current_time - network.last_non_operational_times[current_event.machine])

                    # Update time stamped maintenance cost statistic
                    network.maintenance_costs_stamped[current_event.machine].append((network.maintenance_costs[current_event.machine] + network.total_downtime_costs[current_event.machine], network.current_time))

                    # Updating total_repair_time statistic
                    network.total_repair_time += network.current_time - network.last_repair_time

                # Setting the current machine's state to "FUNCTIONAL" and resetting its degradation level
                network.machines[current_event.machine].state = "FUNCTIONAL"
                network.machines[current_event.machine].degradation = 0

                # If there is dependency update all the arrival intensities
                if network.dependency is True:

                    network.handle_dependency(network.lambdas)

                # Schedule the next degradation event for the machine
                network.schedule_degradation(current_event.machine)

                # Handle the fse action
                network.handle_fse_action_greedy(current_event)

    # Printing the statistics
    print("Correct Tours:", collections.Counter(network.correct_tour))
    print("Wrong Tours:", collections.Counter(network.detour))
    remaining = [list(map(itemgetter(i), network.remaining_times)) for i in range(network.nr_machines)]
    print("Average Remaining Times:", [sum(i) / len(i) for i in remaining])
    print("Number of Repairs:", network.number_of_repairs)
    print("Mean Jumps: ", [sum(i) / len(i) for i in network.jump_sizes])
    print("Mean Interarrival Times:", [sum(i) / len(i) for i in network.interarrivals])
    print("Number of failed machines at decision points:", collections.Counter(network.failed_count))
    print("Idle Percentage of FSE:", (network.total_idle_time / (max_time - network.warm_up_time)) * 100)
    print("Travel Percentage of FSE:", (network.total_travel_time / (max_time - network.warm_up_time)) * 100)
    print("Repair Percentage of FSE:", (network.total_repair_time / (max_time - network.warm_up_time)) * 100)
    print("Percentage Check:", (network.total_idle_time + network.total_travel_time + network.total_repair_time) / (max_time - network.warm_up_time))
    print("REMOVED EVENT COUNT:", network.removed_event_count)
    return network.non_operational_times, network.response_times, network.maintenance_costs, network.total_downtime_costs, \
        network.non_operational_times_stamped, network.maintenance_costs_stamped

# Method to read the input files with a given name
# Returns several input variables for the machines: failure thresholds, arrival intensities,
# alpha and beta shape parameters, maintenance durations, maintenance and downtime costs


def read_input(machine_input_file, travel_time_file):

    # Initialise all input variables to return except travel times
    nr_machines = 0
    list_thresholds = []
    list_lambdas = []
    list_alphas = []
    list_betas = []
    list_preventative_times = []
    list_corrective_times = []
    list_preventative_costs = []
    list_corrective_costs = []
    list_downtime_costs = []

    # Initialise linecount
    linecount = 0

    # Open the machine input file with given filename
    with open(machine_input_file, 'r') as file:

        # Iterate over each line
        for line in file:

            # Get the number of stations from the line length
            if linecount == 0:

                nr_machines = len(line.split())

                # Initialize travel times matrix
                list_travel_times = [[0] * nr_machines for i in range(nr_machines)]

            # Iterate over each word in the line
            for word in line.split():

                # Convert the word to a float
                number = float(word)

                # Line 0 contains the failure thresholds
                if linecount == 0:
                    list_thresholds.append(number)

                # Line 1 contains the lambdas
                elif linecount == 1:
                    list_lambdas.append(number)

                # Line 2 contains the alpha shape parameters
                elif linecount == 2:
                    list_alphas.append(number)

                # Line 3 contains the beta shape parameters
                elif linecount == 3:
                    list_betas.append(number)

                # Line 4 contains the preventative repair times
                elif linecount == 4:
                    list_preventative_times.append(number)

                # Line 5 contains the corrective repair times
                elif linecount == 5:
                    list_corrective_times.append(number)

                # Line 6 contains the preventative repair costs
                elif linecount == 6:
                    list_preventative_costs.append(number)

                # Line 7 contains the corrective repair costs
                elif linecount == 7:
                    list_corrective_costs.append(number)

                # Line 8 contains the downtime costs
                elif linecount == 8:
                    list_downtime_costs.append(number)

            # Iterate linecount
            linecount += 1

    # Initialise wordcount/linecount
    linecount = 0
    wordcount = 0

    # Open the travel time file with given filename
    with open(travel_time_file, 'r') as file:

        # Iterate over each line
        for line in file:

            # Iterate over each word in the line
            for word in line.split():

                # Convert the word to a float
                number = float(word)

                list_travel_times[linecount][wordcount] = number

                # Iterate wordcount
                wordcount += 1

            # Iterate linecount
            wordcount = 0
            linecount += 1

    # Return all collected input variables
    return (nr_machines, list_thresholds, list_lambdas, list_alphas, list_betas, list_preventative_times,
            list_corrective_times, list_preventative_costs, list_corrective_costs, list_downtime_costs,
            list_travel_times)


def calc_conf_interval_list(lst, nr_runs):
    '''
    Method to calculate confidence intervals of result variables
    Note that we use 1.96, thus 95% confidence interval
    '''

    sd = np.std(lst)
    m = np.mean(lst)

    halfwidth = 1.96 * sd / math.sqrt(nr_runs)
    lower = m - halfwidth
    upper = m + halfwidth

    return round(lower, 4), round(upper, 4), round(m, 4)


if __name__ == "__main__":

    # Print the starting time of the simulation
    dateTimeObj = datetime.now()
    print("Current time is: ")
    print(dateTimeObj.hour, ':', dateTimeObj.minute, ':', dateTimeObj.second)

    # Run simulation with given input files
    main("input.txt", "matrix.txt")

    # Print the ending time of the simulation
    dateTimeObj = datetime.now()
    print("Finishing time is: ")
    print(dateTimeObj.hour, ':', dateTimeObj.minute, ':', dateTimeObj.second)
