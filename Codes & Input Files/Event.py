# Event class based on the SwitchingServerQueue example provided by the course.
# Contains the type of event, the time of the event and the machine related to that event:
# END_MAINTENANCE Event => The machine that the FSE has finished repairing
# ARRIVAL Event => The machine that the FSE has arrived
# DEGRADATION Event => The machine that is going to depreciate

class Event:

    def __init__(self, typ, time, machine):

        self.typ = typ  # END_MAINTENANCE, ARRIVAL or DEGRADATION
        self.time = time  # real positive number
        self.machine = machine  # machine number

    def __lt__(self, other):
        return self.time < other.time

    def __str__(self):
        return self.typ + " event " + ' at t = ' + str(self.time)
