#!/usr/bin/python3
from statistics import *

def get_data_from_parallel(filename):
    f = open(filename, "r")
    file_as_string = f.read()
    file_as_array = file_as_string.split("\n")

    shared_allocated_gpu_time_str = file_as_array[0].split(":")[1].strip()
    shared_allocated_gpu_time     = float(shared_allocated_gpu_time_str)

    cpu_time_str = file_as_array[3].split(":")[2].strip()
    cpu_time     = float(cpu_time_str)

    gpu_time_str = file_as_array[4].split(":")[2].strip()
    gpu_time     = float(gpu_time_str)

    return [cpu_time, shared_allocated_gpu_time + gpu_time]

def get_data_from_sequential(filename):
    f = open(filename, "r")
    file_as_string = f.read()
    file_as_array = file_as_string.split("\n")

    cpu_time_str = file_as_array[1].split(":")[2].strip()
    cpu_time     = float(cpu_time_str)

    return cpu_time

def get_all_sequential_data(mesh, n):
    cpu_times = []
    sequential_name = "sequential_" + str(mesh) + "_"
    extention = ".txt"

    for i in range(1,n+1):
        filename = sequential_name + str(i) + extention
        cpu_time = get_data_from_sequential(filename)
        cpu_times.append(cpu_time)

    return cpu_times

def get_all_parallel_data(mesh, n):
    cpu_times = []
    gpu_times = []
    parallel_name = "parallel_" + str(mesh) + "_"
    extention = ".txt"

    for i in range(1,n+1):
        filename = parallel_name + str(i) + extention
        times = get_data_from_parallel(filename)
        cpu_times.append(times[0])
        gpu_times.append(times[1])

    return cpu_times, gpu_times

def main():
    mesh_numbers = (240, 960, 2160, 4000)
    cpu_sequential_times = []
    cpu_parallel_times = []
    gpu_times = []
    n = len(mesh_numbers)
    datas = 30

    for i in range(n):
        mesh = mesh_numbers[i]
        
        times = get_all_sequential_data(mesh, datas)
        cpu_sequential_times.append(times)

        times = get_all_parallel_data(mesh, datas)
        cpu_parallel_times.append(times[0])
        gpu_times.append(times[1])
    
        print("CPU sequencial ", mesh, "média = ", mean(cpu_sequential_times[i]), " mediana = ", median(cpu_sequential_times[i]), "variancia = ", variance(cpu_sequential_times[i]))
        print("CPU parallel   ", mesh, "média = ", mean(cpu_parallel_times[i]), " mediana = ", median(cpu_parallel_times[i]), "variancia = ", variance(cpu_parallel_times[i]))
        print("GPU parallel   ", mesh, "média = ", mean(gpu_times[i]), " mediana = ", median(gpu_times[i]), "variancia = ", variance(gpu_times[i]))

if __name__ == "__main__":
    main()
