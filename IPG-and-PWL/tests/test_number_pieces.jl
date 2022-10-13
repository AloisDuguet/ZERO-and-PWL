using LinA
#alpha = 1 # not working in an expr
expr_1D = :(1*(1/(sqrt(1-x))-1))
err = Absolute(0.05) # sufficiently small?
t1 = 0
t2 = 0.92 # B_i = 2.5 and h_i(0.92) ~= 2.54
pwl = LinA.exactLin(expr_1D,t1,t2,err)
# solution is 25 pieces


#= solution for Absolute(0.05) and t2=0.92 with 4 pieces
4-element Vector{LinA.LinearPiece}:
 0.9722107017380843 x -0.05 from 0 to 0.6021611009198456
 3.431022517811856 x -1.5306008301217076 from 0.6021611009198456 to 0.8116516499227371
 9.68639614420093 x -6.607785154863575 from 0.8116516499227371 to 0.8996206863005184
 23.515841644057236 x -19.04904040659992 from 0.8996206863005184 to 0.92=#
