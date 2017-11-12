#!/usr/bin/python3
from statistics import *
import re

float_regex = "([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)"
keys  = ("SHARED", "GHMATECE", "RIGID", "GHMATECD", "LINSOLVE", "INTEREC1")
mesh_numbers = (240, 960, 2160, 4000, 14400)
modes = ("cpu", "gpu", "gpu_sing")
threads = (1, 8, 24, 48)
executions = 30

def get_data_from_file(filename):
    result = {}
    f = open(filename, "r")
    regexes = (
        "Tempo gasto alocando os vetores compartilhados:\s*" + float_regex,
        "GHMATECE: Tempo na ...:\s*" + float_regex,
        "GHMATECE: Corpo rigido:\s*" + float_regex,
        "GHMATECD: Tempo na ...:\s*" + float_regex,
        "LINSOLVE: Tempo na ...:\s*" + float_regex,
        "INTEREC1: Tempo na ...:\s*" + float_regex
    )

    lines = f.read().split("\n")
    for line in lines:
        for i in range(len(keys)):
            key = keys[i]
            regex = regexes[i]

            match = re.search(regex, line)
            if (match):
                result[key] = float(match.group(1))
                break

    return result


def get_all_data(mode, mesh, threads, n):
    name = "results/results_" + mode + "_" + str(mesh) + "_" + str(threads) + "_"
    extention = ".txt"

    times = {}
    for key in keys:
        times[key] = []

    for i in range(1,n+1):
        filename = name + str(i) + extention
        time = get_data_from_file(filename)
        for key in time:
            times[key].append(time[key])
        
    return times

def generate_r_data(values, subroutine_name):
    r_format1 = "{} <- cbind(size, \"{}\", c("
    r_format2 = "))"

    for mode in modes:
        for thread in threads:
            string = ""
            for mesh in mesh_numbers:
                string += str(values[mode][thread][mesh][subroutine_name]["MEAN"]) + ","
       
            string = string[:-2]
            final_string = r_format1.format(mode + str(thread), mode + str(thread)) + string + r_format2
            print(final_string)

    string = ""
    for mode in modes:
        for thread in threads:
            for mesh in mesh_numbers:
                string += str(values[mode][thread][mesh][subroutine_name]["STDEV"]) + ","
       
    string = string[:-2]
    final_string = r_format1.format("DP", "DP") + string + r_format2
    print(final_string)

def main():
    statistical_values = {}
    for mode in modes:
        statistical_values[mode] = {}
        for thread in threads:
            statistical_values[mode][thread] = {}
            for mesh in mesh_numbers:
                statistical_values[mode][thread][mesh] = {}
                t = get_all_data(mode, mesh, thread, executions)
                for key in t:
                    if (len(t[key]) == 0):
                        continue
                    m  = mean(t[key])
                    sd = stdev(t[key])
                    
                    statistical_values[mode][thread][mesh][key] = {}
                    statistical_values[mode][thread][mesh][key]["MEAN"] = m
                    statistical_values[mode][thread][mesh][key]["STDEV"] = sd


    generate_r_data(statistical_values, "GHMATECD")


if __name__ == "__main__":
    main()

