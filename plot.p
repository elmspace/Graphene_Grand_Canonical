reset
unset key

file = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Grand_Canonical/RESULTS/MOD1.dat"

set pointsize 2

plot file using 6:($11/4.77) w lp lw 2 pt 6,\
file using 6:($10/4.17) w lp lw 2 pt 6,\
file using 6:($9/2.36) w lp lw 2 pt 6

pause(-1)

