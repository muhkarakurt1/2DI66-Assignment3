# FES class based on the SwitchingServerQueue example provided by the course.
import heapq


class FES:

    def __init__(self):
        self.events = []

    def add(self, event):
        heapq.heappush(self.events, event)

    def next(self):
        return heapq.heappop(self.events)

    def is_empty(self):
        return len(self.events) == 0

    def get_length(self):
        return len(self.events)

    def peek(self):
        return self.events[0]

    def __str__(self):

        s = ''
        sorted_events = sorted(self.events)
        for e in sorted_events:
            s += str(e) + '\n'
        return s
