reset
unset key

set pointsize 2

plot "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Grand_Canonical/RESULTS/MOD_main_delfE_vs_delV.dat" using 6:($11/7.25) w lp lw 2 pt 6,\
"/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Grand_Canonical/RESULTS/MOD_main_delfE_vs_delV.dat" using 6:($10/5.275) w lp lw 2 pt 6,\
"/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Grand_Canonical/RESULTS/MOD_main_delfE_vs_delV.dat" using 6:($9/2.425) w lp lw 2 pt 6

pause(-1)

