def list_to_string(l):
    # return a string format of list l
    s = ""
    for i in range(len(l)):
        #-0.0082151802835142, -0.00033003269527398515 pour %f
        #-0.008215286798986199, -0.0003300554362795083
        s += "%f"%(l[i])
        if i != len(l)-1:
            s+= " "
    return s

def save_results_SGM(filename, ne, profits, S, n_iter, cpu_time):
    # write in files the output of a run of Compute_NE.IterativeSG_NOT_DFS

    # open file
    file = open(filename, "a")

    # write NE
    file.write(list_to_string(ne)+"\n")

    # write other informations
    file.write("profits "+list_to_string(profits)+"\n")
    file.write("%i iterations\n"%n_iter)
    file.write("%f seconds\n"%cpu_time)

    # write S
    for p in range(len(S)):
        file.write("strategies of player %i:\n"%(p+1))
        for i in range(len(S[p])):
            file.write(list_to_string(S[p][i])+"\n")


    file.write("\n") # to make space between this solution and the others
    file.close()
    return 0
