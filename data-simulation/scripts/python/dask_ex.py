if __name__ == "__main__":
    from dask.distributed import Client

    client = Client()

    # Define your own code
    def f(x):
        return x + 1

    # Run your code in parallel
    futures = client.map(f, range(100))
    results = client.gather(futures)
    print(results)
    client.close()
