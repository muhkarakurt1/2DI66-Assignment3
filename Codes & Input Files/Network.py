import numpy as np

from FSE import FSE
from FES import FES
from Machine import Machine
from Event import Event
import operator


class Network:

    def __init__(self, nr_machines, policy, list_thresholds, list_lambdas, list_alphas, list_betas,
                 list_preventative_times, list_corrective_times, list_preventative_costs,
                 list_corrective_costs, list_downtime_costs, list_travel_times, dependency, warm_up_time):

        self.nr_machines = nr_machines
        self.policy = policy
        self.thresholds = list_thresholds

        # Both base lambdas and lambda's with dependency are created
        # If there is no dependency lambda's with dependency is equal to base lambdas and don't change during the simulation
        self.lambdas = list_lambdas.copy()
        self.dependency_lambdas = list_lambdas.copy()

        self.alphas = list_alphas
        self.betas = list_betas
        self.preventative_times = list_preventative_times
        self.corrective_times = list_corrective_times
        self.preventative_costs = list_preventative_costs
        self.corrective_costs = list_corrective_costs
        self.downtime_costs = list_downtime_costs
        self.travel_times = list_travel_times
        self.dependency = dependency

        self.machines = []

        # Initializing arrays/lists that will be used for calculating output variables

        # Stores total non-operational time for each machine during the whole simulation
        self.non_operational_times = np.zeros(nr_machines)

        # Stores response times for each machine during the whole simulation
        self.response_times = [[] for i in range(nr_machines)]

        # Stores total maintenance(preventive and corrective) costs for each machine during the whole simulation
        self.maintenance_costs = np.zeros(nr_machines)

        # Stores total downtime costs for each machine during the whole simulation
        self.total_downtime_costs = np.zeros(nr_machines)

        # Initializing time-stamped output variables that are going to be used in plotting
        self.maintenance_costs_stamped = [[] for i in range(nr_machines)]
        self.non_operational_times_stamped = [[] for i in range(nr_machines)]

        # Array to record most recent non operational starting time
        # (Either because of a failure or because of a start of the preventative maintenance)
        self.last_non_operational_times = np.zeros(nr_machines)

        # Initializing additional stats

        # Initializing a list that stores interarrival times of degradation events for each machine
        self.interarrivals = [[] for i in range(nr_machines)]

        # Stats to keep track of wasted or not wasted FSE travels
        self.correct_tour = []
        self.detour = []

        # Statistic to keep track of remaining degradation levels to failure at each decision point
        self.remaining_times = []

        # Statistic to keep track of remaining degradation levels to failure at each decision point
        self.number_of_repairs = np.zeros(nr_machines)

        # Statistic to keep track of jump sizes for each machine
        self.jump_sizes = [[] for i in range(nr_machines)]

        # Statistic to keep track of the last start of idle,travel,repair time for FSE
        self.last_idle_time = 0
        self.last_travel_time = 0
        self.last_repair_time = 0

        # Statistic to keep track of the total idle,travel,repair time for FSE
        self.total_idle_time = 0
        self.total_travel_time = 0
        self.total_repair_time = 0

        # Statistic to keep track of the total number of failed machines at each decision point
        self.failed_count = []

        # Statistic to keep track of the total number of skipped degradation events because of start of a repair on that machine
        self.removed_event_count = 0

        # Time statistics
        self.current_time = 0
        self.warm_up_time = warm_up_time  # Until the warm-up time is reached the output variables are not updated

        # Initializing the field service engineer
        self.fse = FSE()

        # Initializing FES to store all events
        self.fes = FES()

        for machine in range(nr_machines):

            # Adding the machines with their respective attributes to the list
            self.machines.append(Machine(machine, list_thresholds[machine], list_lambdas[machine],
                                         list_alphas[machine], list_betas[machine],
                                         list_preventative_times[machine], list_corrective_times[machine],
                                         list_preventative_costs[machine], list_corrective_costs[machine],
                                         list_downtime_costs[machine]))

            # Scheduling the first degradation for each machine
            self.schedule_degradation(machine)

    def schedule_arrival(self, initial_machine, target_machine):

        '''
        Function to schedule the arrival of field service engineer event from initial machine to target machine

        '''

        # Retrieve the travel time between initial machine and destination machine
        travel_time = self.travel_times[initial_machine][target_machine]

        # Schedule arrival event at current_time + travel_time
        event = Event("ARRIVAL", self.current_time + travel_time, target_machine)

        # Add scheduled event to the FES
        self.fes.add(event)

    def schedule_degradation(self, machine_nr):

        '''
        Function to schedule degradation event for given machine

        '''
        # Compute the time between degradation events for given machine
        # Interarrival times of degradation events are always calculated with lambda's with dependency
        interarrival_time = np.random.exponential(1 / self.dependency_lambdas[machine_nr])

        # Updating the interarrivals statistic
        self.interarrivals[machine_nr].append(interarrival_time)

        # Schedule degradation event at current_time + interarrival_time
        event = Event("DEGRADATION", self.current_time + interarrival_time, machine_nr)

        # Add scheduled event to the FES
        self.fes.add(event)

    def schedule_end_maintenance(self, machine_nr, repair_type):

        '''
        Function to schedule end of maintenance event for the field service engineer

        '''
        if repair_type == "CORRECTIVE":

            # Schedule end of maintenance event at current_time + corresponding corrective maintenance time
            event = Event("END_MAINTENANCE", self.current_time + self.corrective_times[machine_nr], machine_nr)

            # Add scheduled event to the FES
            self.fes.add(event)

        elif repair_type == "PREVENTIVE":

            # Schedule end of maintenance event at current_time + corresponding preventative maintenance time
            event = Event("END_MAINTENANCE", self.current_time + self.preventative_times[machine_nr], machine_nr)

            # Add scheduled event to the FES
            self.fes.add(event)

    def find_max_downtime(self):

        '''
        Function that finds the machine with the maximum downtime cost among failed machines
        (If there are multiple failed machines with the same downtime cost ties are resolved randomly)

        Returns selected machine and total number of failed machines at that decision point
        '''

        machine_number_and_cost = [(machine.machine_nr, machine.down_cost) for machine in self.machines if machine.state == "FAILED"]

        # If there are no failed machines at that decision point return -1 and 0(no failed machine)
        if len(machine_number_and_cost) == 0:

            return -1, 0

        # Calculating the maximum downtime cost among all failed machines
        max_downtime_penalty = max(machine_number_and_cost, key=lambda x: x[1])[1]

        # Selecting among all machines with the maximum downtime cost randomly
        selected_machine = np.random.choice([t[0] for t in machine_number_and_cost if t[1] == max_downtime_penalty])

        # Return the selected machine and the total number of failed machines at that decision point
        return selected_machine, len(machine_number_and_cost)

    def find_biggest_degradation(self):

        '''
        Function that finds the machine with the biggest degradation among non-healthy machines

        Ties are first resolved by travel times and then further ties are resolved randomly)

        Returns selected machine and if that machine is FAILED or not at that decision point
        '''

        current_location = self.fse.location
        machine_number_degradation_traveltime = [(machine.machine_nr, machine.degradation, self.travel_times[current_location][machine.machine_nr]) for machine in self.machines if machine.degradation != 0]

        # If all the machines are healthy at that decision point return -1 and True(Functional Machine)
        if len(machine_number_degradation_traveltime) == 0:

            return -1, True

        # Sorting the list of tuples based on two criteria
        # (Degradation element is descending, travel time element is ascending)
        sorted_list = sorted(machine_number_degradation_traveltime, key=operator.itemgetter(2))
        sorted_list = sorted(sorted_list, key=operator.itemgetter(1), reverse=True)

        # Getting the max_degradation of the machine with max_degradation
        max_degradation = sorted_list[0][1]

        # Getting the min_travel_time of the machine with max_degradation
        min_travel_time = sorted_list[0][2]

        # Selecting among these machines randomly
        selected_machine = np.random.choice([t[0] for t in machine_number_degradation_traveltime if t[1] == max_degradation and t[2] == min_travel_time])

        # Returns selected machine and if that machine is FAILED or not at that decision point
        return selected_machine, self.machines[selected_machine].state == "FUNCTIONAL"

    def handle_fse_action_greedy(self, current_event):

        '''
        Function that handles all the required FSE actions for greedy policy at a decision point

        Takes current event as input parameter
        '''

        # Selecting the machine with the biggest degradation among non-healthy machines
        # and returning the condition of that machine
        selected_machine, functional_or_not = self.find_biggest_degradation()

        # If all the machines are healthy (have 0 degradation) FSE remains idle
        if selected_machine == -1:

            self.fse.state = "IDLE"

            # Update the last idle time
            self.last_idle_time = self.current_time

        # If that selected machine is the FSE's current location start one of the maintenance events
        elif selected_machine == self.fse.location:

            # If FSE starts maintenance at the same machine FSE set out to increase correct_tour statistic
            if current_event.typ == "ARRIVAL" and self.current_time > self.warm_up_time:

                # Update the correct_tour statistic if the simulation time is above warm-up time
                self.correct_tour.append(current_event.machine)

            # If the current state of fse is IDLE update the total_idle_time
            if self.fse.state == "IDLE" and self.current_time > self.warm_up_time:

                self.total_idle_time += self.current_time - self.last_idle_time

            # Set the field service engineer's state to "REPAIR" and update last repair time statistic
            self.fse.state = "REPAIR"
            self.last_repair_time = self.current_time

            # If the selected machine is still functional
            if functional_or_not is True:

                # Update the last non-operational time of the machine
                self.last_non_operational_times[current_event.machine] = self.current_time

                if self.current_time > self.warm_up_time:

                    # Increasing the total maintenance cost for the machine
                    self.maintenance_costs[current_event.machine] += self.preventative_costs[current_event.machine]

                # Setting the degradation level to failure threshold to indicate non-operational machine
                self.machines[selected_machine].degradation = self.machines[selected_machine].failure_threshold
                self.machines[selected_machine].state = "FAILED"

                # If there is dependency update all the arrival intensities
                if self.dependency is True:

                    self.handle_dependency(self.lambdas)

                # Schedule the end of repair event for preventative maintenance
                self.schedule_end_maintenance(current_event.machine, "PREVENTIVE")

            elif functional_or_not is False:

                if self.current_time > self.warm_up_time:

                    # Calculate the response time for the machine and update the response_times statistic
                    response_time = self.current_time - self.last_non_operational_times[current_event.machine]
                    self.response_times[current_event.machine].append(response_time)

                    # Increasing the total maintenance cost for the machine
                    self.maintenance_costs[current_event.machine] += self.corrective_costs[current_event.machine]

                    # Update time stamped maintenance cost statistic
                    self.maintenance_costs_stamped[current_event.machine].append((self.maintenance_costs[current_event.machine] + self.total_downtime_costs[current_event.machine], self.current_time))

                # Schedule the end of repair event for corrective maintenance
                self.schedule_end_maintenance(current_event.machine, "CORRECTIVE")

        # If that selected machine is other than the current machine
        else:

            # If the current event is arrival it means FSE have set out to the current machine to fix it and then decided to
            # travel to some other machine when he arrives there, in this case detour statistic is updated
            if current_event.typ == "ARRIVAL" and self.current_time > self.warm_up_time:

                self.detour.append((self.fse.location, selected_machine))

            # If the current state of fse is IDLE update the total_idle_time
            if self.fse.state == "IDLE" and self.current_time > self.warm_up_time:

                # Updating the total idle time
                self.total_idle_time += self.current_time - self.last_idle_time

            # If that selected machine isn't the current machine start travelling
            self.fse.state = "TRAVEL"

            # Update the last travel time statistic
            self.last_travel_time = self.current_time

            # Schedule the arrival event for the field service engineer
            self.schedule_arrival(self.fse.location, selected_machine)

            # Updating the remaining_times statistic
            self.remaining_times.append([10 - machine.degradation for machine in self.machines])

    def handle_fse_action_reactive(self, current_event):

        '''
        Function that handles all the required FSE actions for reactive policy at a decision point

        Takes current event as input parameter
        '''

        # Selecting the machine with the maximum downtime cost among failed machines
        selected_machine, failed_machines_count = self.find_max_downtime()

        # Update the failed_count statistic if the simulation time is above warm-up time
        if self.current_time > self.warm_up_time:

            # Update the failed machine count statistic
            self.failed_count.append(failed_machines_count)

        # If no machine is failed at the moment FSE remains idle
        if selected_machine == -1:

            self.fse.state = "IDLE"

            # Update the last idle time statistic
            self.last_idle_time = self.current_time

        # If that selected machine is the current machine start the corrective maintenance
        elif selected_machine == self.fse.location:

            # If FSE starts corrective maintenance at the same machine FSE set out to increase correct_tour statistic
            if current_event.typ == "ARRIVAL" and self.current_time > self.warm_up_time:

                # Update the correct_tour statistic if the simulation time is above warm-up time
                self.correct_tour.append(current_event.machine)

            # If the current state of fse is IDLE update the total_idle_time
            if self.fse.state == "IDLE" and self.current_time > self.warm_up_time:

                # Updating the total idle time
                self.total_idle_time += self.current_time - self.last_idle_time

            # Set the field service engineer's state to "REPAIR"
            self.fse.state = "REPAIR"
            self.last_repair_time = self.current_time

            # Update the performance measures if the simulation time is above warm-up time
            if self.current_time > self.warm_up_time:

                # Calculate the response time for the machine and update the response_times statistic
                response_time = self.current_time - self.last_non_operational_times[current_event.machine]
                self.response_times[current_event.machine].append(response_time)

                # Increasing the total maintenance cost for the machine
                self.maintenance_costs[current_event.machine] += self.corrective_costs[current_event.machine]

                # Update time stamped maintenance cost statistic
                self.maintenance_costs_stamped[current_event.machine].append((self.maintenance_costs[current_event.machine] + self.total_downtime_costs[current_event.machine], self.current_time))

            # Schedule the end of repair event
            self.schedule_end_maintenance(current_event.machine, "CORRECTIVE")

        # If that selected machine is other than the current machine
        else:

            # If the current event is arrival it means FSE have set out to the current machine to fix it and then decided to
            # travel to some other machine when he arrives there, in this case detour statistic is updated
            if current_event.typ == "ARRIVAL" and self.current_time > self.warm_up_time:

                self.detour.append((self.fse.location, selected_machine))

            # If the current state of fse is IDLE update the total_idle_time
            if self.fse.state == "IDLE" and self.current_time > self.warm_up_time:

                # Updating the total idle time
                self.total_idle_time += self.current_time - self.last_idle_time

            # Set the state of field service engineer to "TRAVEL"
            self.fse.state = "TRAVEL"
            self.last_travel_time = self.current_time

            # Schedule the arrival event for the field service engineer
            self.schedule_arrival(self.fse.location, selected_machine)

            # Updating the remaining_times statistic
            self.remaining_times.append([10 - machine.degradation for machine in self.machines])

    def handle_degradation(self, current_event):

        '''
        Function that increases current event's machine's degradation level

        Takes current event as input parameter
        '''

        # Retrieving the alpha and beta parameters for the current machine
        alpha = self.alphas[current_event.machine]
        beta = self.betas[current_event.machine]

        # Finding the jump size(sampling from beta distribution)
        jump_size = np.random.beta(alpha, beta)

        # Updating the jump_sizes statistic
        self.jump_sizes[current_event.machine].append(jump_size)

        # Updating the degradation level for the current machine
        self.machines[current_event.machine].degradation += jump_size

    def handle_dependency(self, base_lambdas):

        '''
        Function that updates all the dependency lambdas for the network
        (When a machine fails or a preventive maintenance is started)

        Takes base_lambdas as input parameter
        '''

        # Find operational and non-operational machines in the network
        non_operational_machines = [machine.machine_nr for machine in self.machines if machine.degradation == machine.failure_threshold]
        operational_machines = [machine.machine_nr for machine in self.machines if machine.degradation < machine.failure_threshold]

        # For non operational machines set lambda to base lambda
        for non_operational_machine in non_operational_machines:

            self.dependency_lambdas[non_operational_machine] = 0

        # If there is no operational machine do nothing
        if len(operational_machines) == 0:

            pass

        # If there are operational machines add non_operational machines' load to them
        elif len(operational_machines) != 0:

            for operational_machine in operational_machines:

                # Calculating additional arrival intensity for each machine that is going to come from non operational machines

                additional_arrival_intensity = (1 / len(operational_machines)) * \
                    sum([base_lambdas[non_operational_machine] for non_operational_machine in non_operational_machines])

                # Adding calculated additional arrival intensity to each machines base lambda
                self.dependency_lambdas[operational_machine] = base_lambdas[operational_machine] + additional_arrival_intensity
