# Machine class that stores given machine parameters, real time degradation level and the state of the machine


class Machine:

    def __init__(self, machine_nr, failure_threshold, arrival_intensity, alpha, beta,
                 prev_time, correct_time, prev_cost, correct_cost, down_cost):

        self.machine_nr = machine_nr
        self.failure_threshold = failure_threshold
        self.arrival_intensity = arrival_intensity
        self.alpha = alpha
        self.beta = beta
        self.prev_time = prev_time
        self.correct_time = correct_time
        self.prev_cost = prev_cost
        self.correct_cost = correct_cost
        self.down_cost = down_cost

        self.degradation = 0
        self.state = "FUNCTIONAL"  # FUNCTIONAL or FAILED
