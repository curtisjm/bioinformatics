import ray
import pandas as pd

@ray.remote
def modify():
    ray.get(ray_df).at[0, 'Name'] = 'joe'

@ray.remote
def printer():
    print(ray.get(ray_df))

ray.init()
print(int(ray.available_resources()["CPU"]))

data = [['tom', 10], ['nick', 15], ['juli', 14]]
df = pd.DataFrame(data, columns=['Name', 'Age'])

ray_df = ray.put(df)

# printer.remote()
# modify.remote()
# printer.remote()
