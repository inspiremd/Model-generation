'''
date: 10/26/2019

taken from
sources: https://gist.github.com/fspaolo/51eaf5a20d6d418bd4d0
'''

from impress_md import interface_functions

import os
import sys
import numpy as np
from mpi4py import MPI
from queue import Queue

WORKTAG = 1
DIETAG = 0


class Work(object):
    def __init__(self):
        # importat: sort by file size in decreasing order!
        q = Queue()
        with open('file_list.txt', 'r') as f:
            paths_to_run_mmgbsa = map(lambda x: x.strip(), f.readlines())
        for f in paths_to_run_mmgbsa:
            q.put(f)
        self.work = q

    def get_next(self):
        if self.work.empty():
            return None
        return self.work.get()


def do_work(work):
    print("running", work)
    interface_functions.RunMinimization(work, work, True)
    print("done with ", work)

# def process_result(result):
#     pass


def master(comm):
    num_procs = comm.Get_size()
    status = MPI.Status()

    # generate work queue
    wq = Work()

    # Seed the slaves, send one unit of work to each slave (rank)
    for rank in range(1, num_procs):
        work = wq.get_next()
        comm.send(work, dest=rank, tag=WORKTAG)

    # Loop over getting new work requests until there is no more work to be done
    while True:
        work = wq.get_next()
        if not work: break

        # Receive results from a slave
        result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        # process_result(result)

        # Send the slave a new work unit
        comm.send(work, dest=status.Get_source(), tag=WORKTAG)

    # No more work to be done, receive all outstanding results from slaves
    for rank in range(1, num_procs):
        result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        # process_result(result)

    # Tell all the slaves to exit by sending an empty message with DIETAG
    for rank in range(1, num_procs):
        comm.send(0, dest=rank, tag=DIETAG)


def slave(comm):
    my_rank = comm.Get_rank()
    status = MPI.Status()

    while True:
        # Receive a message from the master
        work = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)

        # Check the tag of the received message
        if status.Get_tag() == DIETAG: break

        # Do the work
        do_work(work)

    # Send the result back
    comm.send(1, dest=0, tag=0)


def main():
    comm = MPI.COMM_WORLD
    my_rank = comm.Get_rank()
    my_name = MPI.Get_processor_name()
    # comm.Barrier()
    # start = MPI.Wtime()

    if my_rank == 0:
        master(comm)
    else:
        slave(comm)

    # comm.Barrier()
    # end = MPI.Wtime()
    # print 'time:', end - start


if __name__ == '__main__':
    main()