import os

script_dir = os.path.dirname(os.path.abspath(__file__))
top_dir = os.path.abspath(script_dir + "/..")
data_dir = top_dir + "/data"
paper_dir = top_dir + "/hipc16"
cmd_exe = script_dir + "/svd"

num_instances = 2
tol = 1e-15
lamda = 0.0
topk = 4

sym_mat_sizes = [500, 750, 1000, 1250, 1500, 1750, 2000]
mat_sizes = [[500, 300], [750, 500], [1000, 700], [1250, 900], [1500, 1300], [1750, 1500], [2000,1800]]
solvers = [1,2,3,4,5,6,7,8,9,10]

sym_mat_sizes = [30, 50]
mat_sizes = [[30, 25], [50, 20]]
solvers = [9]


