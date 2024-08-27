from ray.util import ActorPool
import ray
import pandas as pd
import numpy as np

@ray.remote
class GlobalStateActor():
    def __init__(self):
        data = [['tom', 10], ['nick', 15], ['juli', 14]]
        self.df = pd.DataFrame(data, columns=['Name', 'Age'])
        self.rng = np.random.default_rng()
    
    def get_df(self):
        return self.df

    def get_rng(self):
        return self.rng.random()
    
    def mutate_df(self, row, col, value):
        self.df.at[row, col] = value

@ray.remote
class SimulationActor():
    def __init__(self, global_state):
        self.global_state = global_state
    
    def modify(self, row, col, value):
        self.global_state.mutate_df.remote(row, col, value)
    
    def printer(self):
        print(ray.get(self.global_state.get_df.remote()))
    
    def random(self):
        return self.global_state.get_rng.remote()
    
    def p(self, i):
        print(i)


    
ray.init()

gs = GlobalStateActor.remote()
sa1 = SimulationActor.remote(gs)
sa2 = SimulationActor.remote(gs)

pool = ActorPool([sa1, sa2])
list(pool.map(lambda actor, i: actor.p.remote(i), range(10))) 

# futures = [sa1.modify.remote(0, "Name", "JOE"), sa2.modify.remote(1, "Name", "TOBY"), sa2.printer.remote()]
# ray.get(futures)
# print(ray.get(sa1.random.remote()))

