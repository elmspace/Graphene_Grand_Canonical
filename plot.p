reset
unset key

set pointsize 2

plot "/Users/ashkandehghan/Desktop/Graphene/Computing/GrandCanonical/main_grand_canonical/RESULTS/p_ave.dat" using 3:1 w lp lw 2 pt 6

pause(-1)


   reset
set pm3d
set iso 100
set samp 100
set palette model RGB
set dgrid3d 20,20,1
set pm3d flush begin ftriangles scansforward interpolate 10,5


   unset key
   unset sur
   set hidden3d
   set view map 
   set autoscale
#set size square
   splot "/Users/ashkandehghan/Desktop/Graphene/Computing/main1/omega_test.dat" using 1:2:4
