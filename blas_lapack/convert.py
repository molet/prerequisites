#!/usr/bin/python

import sys
import numpy as np

for arg in sys.argv:
    with open(arg) as f:
        con = f.readlines()

for i in range(np.size(con)):
    # comment
    if con[i][0] == "*":
        con[i]=con[i].replace("*","!",1)
    elif con[i].find("SUBROUTINE") != -1:
        cls="SUBROUTINE"
    elif con[i].find("FUNCTION") != -1:
        cls="FUNCTION"
    elif con[i].strip() == "END":
        con[i]=con[i][:len(con[i])-1]+" "+cls+con[i][len(con[i])-1:]
    # cont
    elif len(con[i]) > 6:
        if con[i][5] == "$":
            con[i]=con[i].replace("$"," ",1)
            con[i-1]=con[i-1][:len(con[i-1])-1]+" &"+con[i-1][len(con[i-1])-1:]
        elif con[i][5] == "+":
            con[i]=con[i].replace("+"," ",1)
            con[i-1]=con[i-1][:len(con[i-1])-1]+" &"+con[i-1][len(con[i-1])-1:]
        elif con[i][5] == "&":
            con[i]=con[i].replace("&"," ",1)
            con[i-1]=con[i-1][:len(con[i-1])-1]+" &"+con[i-1][len(con[i-1])-1:]

i=0
while i < np.size(con):
    if (con[i].find("Intrinsic Functions") != -1) or \
        (con[i].find("External Function") != -1) or \
        (con[i].find("External Subroutines") != -1) or \
        (con[i].find("     from BLAS") != -1) or \
        (con[i].find("     from LAPACK") != -1):
        i+=1
        while con[i].find("..") == -1:
            con[i]=con[i].replace(" ","!",1)
            i+=1
    else:
        i+=1

for i in range(np.size(con)):
    # dorbdb.f
    if con[i].find(". GT.") != -1:
       con[i]=con[i].replace(". GT.",".GT.")
    # dtgex2.f
    if con[i].find("2.0D+01") != -1:
        con[i]=con[i].replace("2.0D+01","20.0_dp")
    if con[i].find("D0") != -1:
        con[i]=con[i].replace("D0","_dp")
    if con[i].find("d0") != -1:
        con[i]=con[i].replace("d0","_dp")
    # dgeqrt2.f, dgeqrt3.f
    if con[i].find("D+00") != -1:
        con[i]=con[i].replace("D+00","_dp")
    if con[i].find("D+0") != -1:
        con[i]=con[i].replace("D+0","_dp")
    if con[i].find("DABS(") != -1:
        con[i]=con[i].replace("DABS(","ABS(")
    if con[i].find("DSQRT(") != -1:
        con[i]=con[i].replace("DSQRT(","SQRT(")
    if con[i].find("DSIGN(") != -1:
        con[i]=con[i].replace("DSIGN(","SIGN(")
    if con[i].find("DLOG(") != -1:
        con[i]=con[i].replace("DLOG(","LOG(")
    if con[i].find("IDNINT(") != -1:
	con[i]=con[i].replace("IDNINT(","INT(")
    if con[i].find("DMIN1(") != -1:
        con[i]=con[i].replace("DMIN1(","MIN(")
    # dlasq2.f, dlasq3.f, dlasq4.f, dlasq5.f, dlasq6.f
    if con[i].find("DMIN1") != -1:
        con[i]=con[i].replace("DMIN1","DMIN_1")
    if con[i].find("DMAX1") != -1:
        con[i]=con[i].replace("DMAX1(","MAX(")
    while con[i].find("DBLE(") != -1:
        p=con[i].find("DBLE(")
        con[i]=con[i].replace("DBLE(","REAL(",1)
        n=1
        for j in range(p+5,len(con[i])):
            if con[i][j:j+1] == "(":
                n+=1
            elif con[i][j:j+1] == ")":
                n-=1
            if n==0:
                con[i]=con[i][:j]+", dp"+con[i][j:]
                break
    if con[i].find("DOUBLE PRECISION") != -1:
        con[i]=con[i].replace("DOUBLE PRECISION","real(dp)")

for i in range(np.size(con)):
    print(con[i]),
