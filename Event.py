# Event class based on the SwitchingServerQueue example provided by the course.
# Contains the type of event, the time of the event and the queue for which the event is.

class Event:
    
    def __init__(self, typ, time, machine):
        
        self.typ = typ  # 
        self.time = time  # real positive number
        self.machine = machine  # machine number
        
    def __lt__(self, other):
        return self.time < other.time

    def __str__(self):
        return self.typ + " of new customer " + ' at t = ' + str(self.time)


