'''
author: Austin Clyde
date: 10/26/2019

sources:
'''


from charm4py import charm

def square(x):
    return x**2

def main(args):
    result = charm.pool.map(square, range(10), chunksize=2)
    print(result)  # prints [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
    exit()

charm.start(main)