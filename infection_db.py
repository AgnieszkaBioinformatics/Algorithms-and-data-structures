from datetime import datetime
from collections import defaultdict
from heapq import heapify, heappop, heappush

class InfectionDB:
    def __init__(self, start_day, people_names, connections, dmax):
        self.start_day = datetime.strptime(start_day, '%Y-%m-%d')
        self.people_names = people_names
        self.connections = connections
        self.dmax = dmax
        self.graph = defaultdict(list)

        for (u, v), w in connections.items():   # O(|V| + |E|)
            self.graph[u].append((v, w))
            self.graph[v].append((u, w))

        self.infections = defaultdict(list)     # {infection_id: [name, date]}

        self.outbreaks = defaultdict(list)      # {first infection: [infections]
        self.largest_outbreak_id = None
        self.infection_id = 0                   # used in add_infection function
        self.number_outbreaks = 0               # used in number_of_outbreaks function
        self.direct_infections = dict()         # used in number_of_direct_infections function
        self.infection_names = set()            # used in add_infection function
    def number_of_outbreaks(self) -> int:       # O(1)
        # returns number of all outbreaks in the db
        return self.number_outbreaks

    def number_of_direct_infections(self, infection_id: int) -> int:    # O(1)
        # returns number of direct infections infected by infection_id
        name = self.infections[infection_id][0]
        if name not in self.direct_infections:
            return 0
        return self.direct_infections[name]

    def number_of_outbreak_infections(self, outbreak_id: int) -> int:   # O(n)
        # returns number of infections in the outbreak identified by outbreak_id
        return len(self.outbreaks[outbreak_id])

    def outbreak_id_of_infection(self, infection_id: int) -> int:   # O(n)
        # returns outbreak_id of specified infection_id
        for key, value in self.outbreaks.items():
            if infection_id in value: return key


    def largest_outbreak_id(self) -> int:    # O(1)
        # returns outbreak_id of the largest outbreak
        return self.largest_outbreak_id

    def outbreak_is_active(self, outbreak_id: int, when: str) -> bool:  # O(n)
        # determine if new infections in the outbreak might occur
        date = datetime.strptime(when, '%Y-%m-%d')
        latest_infection_date = max([self.infections[infection_id][1] for infection_id in self.outbreaks[outbreak_id]])
        return (date - latest_infection_date).days <= self.dmax

    def cleaning(self, who_name: str) -> int:  # O(n)
        # creating new dictionary without previous infections of specified individual
        # so, we'll keep only the last infection
        self.infections = {key: value for key, value in self.infections.items() if who_name not in value}
        return 0

    def add_infection(self, who_name: str, when: str, source_id = None) -> int: # O(n)
        # add new infection to the db
        if who_name not in self.infection_names:  # O(1) according to this stack (for set()):
            self.infection_names.add(who_name)    # https://stackoverflow.com/questions/13884177/complexity-of-in-operator-in-python
        else:
            self.cleaning(who_name)

        when_date = datetime.strptime(when, '%Y-%m-%d')
        infection_id = self.infection_id
        self.infections[infection_id] = [who_name, when_date]
        self.infection_id += 1

        if source_id is None:
            # outbreak_id is the first occurrence of the infection
            outbreak_id = infection_id
            self.outbreaks[infection_id].append(infection_id)
            self.number_outbreaks += 1
        else:
            outbreak_id = self.outbreak_id_of_infection(source_id)
            self.outbreaks[outbreak_id].append(infection_id)

            # if source was specified by the id number
            if isinstance(source_id, int):
                source_name = self.infections[source_id][0]
            else:
                # if source was specified by the name
                source_name = source_id

            if source_name not in self.direct_infections:       # O(n)
                self.direct_infections[source_name] = 1
            else:
                self.direct_infections[source_name] += 1

        # checking if the current outbreak is bigger than the previous one
        if (self.largest_outbreak_id is None or
                len(self.outbreaks[outbreak_id]) > len(self.outbreaks[self.largest_outbreak_id])):
            self.largest_outbreak_id = outbreak_id

        return infection_id

    def dijkstra(self, who_name):      # most probable: O(v + e) x V *
        # modified dijkstra algorithm where we don't analyze the whole graph
        # but only until the distance is not larger than dmax
        distances = dict()
        distances2 = {nid: float("inf") for nid in self.graph} # O(V)
        distances[who_name] = 0
        distances2[who_name] = 0
        visited = set()
        pq = [(who_name, 0)]
        heapify(pq)
        while pq:   # in this implementation it won't be full O(V + E) but the nodes which dist < dmax
            current, length = heappop(pq)
            if current in visited: continue
            visited.add(current)
            for neighbor, weight in self.graph[current]:
                last = (next(reversed(distances)), distances[next(reversed(distances))])[1] # O(1) (in the best case) **
                distance = length + weight
                if last + distance > self.dmax: break
                if distance < distances2[neighbor]:
                    distances[neighbor] = distance
                    distances2[neighbor] = distance
                    heappush(pq, (neighbor, distance))
        return distances
# * v and e are fractions of V and E
# ** according to this stack: https://stackoverflow.com/questions/65540349/time-complexity-of-reversed-in-python-3


    def potential_infection_sources(self, who_name, when):      # O(v + e) x V x n x m
        # return potential sources of infection inferred from the distance on the graph and
        # time duration since the illness started
        sus = self.dijkstra(who_name)       # O((V+E)logV)
        when = datetime.strptime(when, '%Y-%m-%d')
        result = []
        for suspect, dist in sus.items():   # O(m)
            for node, date in self.infections.values():     # O(n)
                if node == suspect:
                    days = (when - date).days
                    if dist < days <= self.dmax:
                        result.append(suspect)
        return result

    def add_infection_inferred_source(self, who_name, when):        # O(v + e) x V  x n^2 x m^2
        # add infection from the inferred source
        # the source was inferred as the last registered probable infection
        sus = self.potential_infection_sources(who_name, when)      # O(v + e) x V x n x m
        last = None
        suspect = 0
        while last != suspect:      # in the best case o(1), in the worst case o(n)
            for suspect in sus:     # O(m)
                last = (next(reversed(self.infections)), self.infections[next(reversed(self.infections))])[1]   #O(1)
                if suspect == last[0]:
                    infection_id = self.add_infection(who_name, when, source_id=suspect)    # O(n)
                    return infection_id


