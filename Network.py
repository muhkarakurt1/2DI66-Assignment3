import numpy as np

from FSE import FSE
from FES import FES
from Machine import Machine
from Event import Event
import operator

class Network:

    def __init__(self, nr_machines, policy, list_thresholds, list_lambdas, list_alphas, list_betas,
                 list_preventative_times, list_corrective_times, list_preventative_costs,
                 list_corrective_costs, list_downtime_costs, list_travel_times, warm_up_time):
        
        self.nr_machines = nr_machines
        self.policy = policy
        self.thresholds = list_thresholds
        self.lambdas = list_lambdas
        self.alphas = list_alphas
        self.betas = list_betas
        self.preventative_times = list_preventative_times
        self.corrective_times = list_corrective_times
        self.preventative_costs = list_preventative_costs
        self.corrective_costs = list_corrective_costs
        self.downtime_costs = list_downtime_costs
        self.travel_times = list_travel_times
        self.machines = []  
        self.interarrivals = [[] for i in range(nr_machines)]
        
        # Initializing arrays/lists that will be used for calculating output variables
        self.non_operational_time = np.zeros(nr_machines)
        self.response_times = [[] for i in range(nr_machines)]
        self.maintenance_costs = np.zeros(nr_machines)
        self.total_downtime_costs = np.zeros(nr_machines)
        
        # Array to record most recent non operational starting time 
        # (either because of a failure or because of a start of the preventative maintenance)
        self.last_non_operational_times = np.zeros(nr_machines)
    
        # Additional stats
        self.correct_tour = []
        self.detour = []
        self.remaining_times = []
        self.number_of_repairs = np.zeros(nr_machines)
        self.jump_sizes = [[] for i in range(nr_machines)]
        self.last_idle_time = 0
        self.total_idle_time = 0
        self.failed_count = []

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
                                         list_preventative_times[machine], list_corrective_costs[machine],
                                         list_preventative_costs[machine], list_corrective_costs[machine],
                                         list_downtime_costs[machine]))

            # Scheduling the first degradation for each machine
            self.schedule_degradation(machine)
    
    # Function to schedule the arrival of field service engineer event from initial machine to destination machine
    def schedule_arrival(self, initial_machine, target_machine):
            
        # Retrieve the travel time between initial machine and destination machine
        travel_time = self.travel_times[initial_machine][target_machine]
    
        # Schedule arrival event at current_time + travel_time
        event = Event("ARRIVAL", self.current_time + travel_time, target_machine)
        
        # Add scheduled event to the FES
        self.fes.add(event)
        
    # Function to schedule degradation event for given machine
    def schedule_degradation(self, machine_nr):
        
        # Compute the time between degradation events for given machine
        interarrival_time = np.random.exponential(1 / self.lambdas[machine_nr])
        
        self.interarrivals[machine_nr].append(interarrival_time)
        # Schedule degradation event at current_time + interarrival_time
        event = Event("DEGRADATION", self.current_time + interarrival_time, machine_nr)
        
        # Add scheduled event to the FES
        self.fes.add(event)
    
    # Function to schedule end of maintenance event for the field service engineer
    def schedule_end_maintenance(self, machine_nr, repair_type):
        
        if repair_type == "CORRECTIVE":
            
            # Schedule end of maintenance event at current_time + corresponding corrective maintenance time
            event = Event("END_MAINTENANCE", self.current_time + self.corrective_times[machine_nr], machine_nr)
        
            # Add scheduled event to the FES
            self.fes.add(event)
            
        elif repair_type == "PREVENTIVE":
                        
            # Schedule end of maintenance event at current_time + corresponding preventative maintenance time
            event = Event("END_MAINTENANCE", self.current_time + self.preventative_costs[machine_nr], machine_nr)
        
            # Add scheduled event to the FES
            self.fes.add(event)
            
    # Function that finds the machine with the maximum downtime penalty among failed machines
    # (ties are resolved randomly)
    def find_max_downtime(self):
        
        machine_number_and_cost = [(machine.machine_nr, machine.down_cost) for machine in self.machines if machine.state == "FAILED"]
        
        if len(machine_number_and_cost) == 0:
            return -1, 0
        
        # Calculating the maximum downtime penalty among all failed machines
        max_downtime_penalty = max(machine_number_and_cost, key=lambda x:x[1])[1]  
        
        # Selecting among all machines with the maximum downtime cost randomly
        selected_machine = np.random.choice([t[0] for t in machine_number_and_cost if t[1] == max_downtime_penalty])
        
        return selected_machine, len(machine_number_and_cost)
        
    # Function that finds the machine with the biggest degradation among non-healthy machines
    # and returns if that machine is FAILED or not (ties are first resolved by travel times and then randomly)  
    def find_biggest_degradation(self):
        
        current_location = self.fse.location
        machine_number_degradation_traveltime = [(machine.machine_nr, machine.degradation, self.travel_times[current_location][machine.machine_nr]) for machine in self.machines if machine.degradation != 0]
        
        # Sorting the list of tuples based on two criteria
        # (degradation element is descending, travel time element is ascending)
        sorted_list = sorted(machine_number_degradation_traveltime, key = operator.itemgetter(2))
        sorted_list = sorted(sorted_list, key = operator.itemgetter(1), reverse=True)
        
        # Getting the max_degradation min_travel_time of the machine with max_degradation
        max_degradation = sorted_list[0][1]
        
        # Getting the min_travel_time of the machine with max_degradation
        min_travel_time = sorted_list[0][2]
        
        # Selecting among these machines randomly
        selected_machine = np.random.choice([t[0] for t in machine_number_degradation_traveltime if t[1] == max_degradation and t[2] == min_travel_time])

        return selected_machine, self.machines[selected_machine].state == "FUNCTIONAL"
    
    def handle_fse_action_greedy(self, current_event):
        
        # Selecting the machine with the biggest degradation among non-healthy machines 
        # and returning the condition of that machine
        selected_machine, functional_or_not = self.find_biggest_degradation()
        
        # If that selected machine is the current machine start one of the maintenance events
        if selected_machine == current_event.machine:
            
            # Updating the total idle time
            self.total_idle_time += self.current_time - self.last_idle_time
                
            #Set the field service engineer's state to "REPAIR"
            self.fse.state = "REPAIR"
    
            # If the selected machine is still functional
            if functional_or_not == True:
                
                self.last_non_operational_times[current_event.machine] = self.current_time
    
                # Increasing the total maintenance cost for the machine
                self.maintenance_costs[current_event.machine] += self.preventative_costs[current_event.machine]
                
                # Schedule the end of repair event for preventative maintenance
                self.schedule_end_maintenance(current_event.machine, "PREVENTIVE")
            
            elif functional_or_not == False:
                
                # Calculate the response time for the machine and update the response_times statistic
                response_time = self.current_time - self.last_non_operational_times[current_event.machine]
                self.response_times[current_event.machine].append(response_time)
    
                # Increasing the total maintenance cost for the machine
                self.maintenance_costs[current_event.machine] += self.corrective_costs[current_event.machine]
                
                # Schedule the end of repair event for corrective maintenance
                self.schedule_end_maintenance(current_event.machine, "CORRECTIVE")
            
        else:
            
            # If that selected machine isn't the current machine start travelling
            self.fse.state = "TRAVEL"
            
            # Schedule the arrival event for the field service engineer
            self.schedule_arrival(current_event.machine, selected_machine)
            
            # Updating the remaining_times statistic
            self.remaining_times.append([10 - machine.degradation for machine in self.machines])
    
    def handle_fse_action_reactive(self, current_event):
        
        # Selecting the machine with the maximum downtime cost among failed machines
        selected_machine, failed_machines_count = self.find_max_downtime()
        
        # Updating the failed_count statistic
        self.failed_count.append(failed_machines_count)
        
        # If that selected machine is the current machine start the corrective maintenance
        if selected_machine == current_event.machine:
            
            if current_event.typ == "ARRIVAL":
                
                self.correct_tour.append(current_event.machine)

            # If the current state of fse is IDLE update the total_idle_time
            if self.fse.state == "IDLE":
                
                # Updating the total idle time
                self.total_idle_time += self.current_time - self.last_idle_time
            
            #Set the field service engineer's state to "REPAIR"
            self.fse.state = "REPAIR"
            
            # Calculate the response time for the machine and update the response_times statistic
            response_time = self.current_time - self.non_operational_time[current_event.machine]
            self.response_times[current_event.machine].append(response_time)

            # Increasing the total maintenance cost for the machine
            self.maintenance_costs[current_event.machine] += self.corrective_costs[current_event.machine]
            
            # Schedule the end of repair event
            self.schedule_end_maintenance(current_event.machine, "CORRECTIVE")
        
        # If no machine is failed at the moment
        elif selected_machine == -1 :
            
            # Setting field service engineer's state to "IDLE"
            self.fse.state = "IDLE"
            self.last_idle_time = self.current_time
            
        # If that selected machine is other than tge current machine
        else:
            
            if current_event.typ == "ARRIVAL":
                
                self.detour.append((current_event.machine, selected_machine))
                
            # Set the state of field service engineer to "TRAVEL"
            self.fse.state = "TRAVEL"
            
            # Schedule the arrival event for the field service engineer
            self.schedule_arrival(current_event.machine, selected_machine)
            
            # Updating the remaining_times statistic
            self.remaining_times.append([10 - machine.degradation for machine in self.machines])

    def handle_degradation(self, current_event):
        
        # Retrieving the alpha and beta parameters for the current machine
        alpha = self.alphas[current_event.machine]
        beta = self.betas[current_event.machine]
        
        # Finding the jump size(sampling from beta distribution)
        jump_size = np.random.beta(alpha, beta)
        
        # Updating the jump_sizes statistic
        self.jump_sizes[current_event.machine].append(jump_size)
        
        # Updating the degradation level for the current machine
        self.machines[current_event.machine].degradation += jump_size

