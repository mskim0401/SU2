import read_su2 as rs

print("Reading flow.dat...")
zones, dfs, variables = rs.read_su2('flow.dat',verbose=1)
list(dfs)