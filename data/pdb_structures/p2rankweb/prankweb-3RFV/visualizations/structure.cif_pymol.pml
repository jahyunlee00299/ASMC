from pymol import cmd,stored

set depth_cue, 1
set fog_start, 0.4

set_color b_col, [36,36,85]
set_color t_col, [10,10,10]
set bg_rgb_bottom, b_col
set bg_rgb_top, t_col      
set bg_gradient

set  spec_power  =  200
set  spec_refl   =  0

load "data/structure.cif", protein
create ligands, protein and organic
select xlig, protein and organic
delete xlig

hide everything, all

color white, elem c
color bluewhite, protein
#show_as cartoon, protein
show surface, protein
#set transparency, 0.15

show sticks, ligands
set stick_color, magenta




# SAS points

load "data/structure.cif_points.pdb.gz", points
hide nonbonded, points
show nb_spheres, points
set sphere_scale, 0.2, points
cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)


stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
lastSTP=stored.list[-1] # get the index of the last residue
hide lines, resn STP

cmd.select("rest", "resn STP and resi 0")

for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))



set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [4611,5314,4186,4613,5099,4614,4615,4599,4600,4601,4602,4605,5125,5126,4609,4703,4671,4697,5319,5321,5324,5326,5408,5406,5329,6043,4899,4900,4908,4916,4918,5398,5325,4621,6065,5377,6045,6000,6048,6053,6064,4887,4888,4892,5315,5312,4336,4170,4174,4175,4179,4166,4167,4598,4337,4351,4338,4342,4345,4449,4458,4459,4149,4589,4595,4152,4153,4327,4328,4330,4331,4590,4593,4155,4183,4185,4525,4466,4460,4465,4467] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.416,0.278,0.702]
select surf_pocket2, protein and id [1038,511,512,513,523,525,526,527,533,615,521,804,1037,812,820,1010,1011,1223,1231,828,829,830,1977,1976,1233,1237,1288,1289,1318,1320,1310,1965,1974,1955,1912,614,370,609,799,800,505,98,1224,1226,1227,95,499,501,502,99,65,247,248,249,67,78,378,377,379,376,507,510,242,254,257,263,243,250,361,372,79,80,91,82,86,89,97,1238,1241,1307,1308,1306] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.902,0.361,0.878]
select surf_pocket3, protein and id [2545,2570,2551,2571,2554,2555,2556,2557,2292,2293,2109,2126,2135,2111,2123,2139,2546,2405,2416,2421,2422,2287,2291,2294,2298,2481,2423,2105,2108,2565,2627,2658,2659,2414,2301] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.380]
select surf_pocket4, protein and id [2864,2856,2872,3054,3055,2569,2577,2133,3354,2873,2874,4018,4020,4021,3352,4001,4004,3999,3956,4009,3294,3285,3267,3268,3270,3275,3277,3280,3281,3282,3333,3363,2844,3271,2142,2546,2135,2111,2139,2141,2545,2570,2561,2843,3081,3082] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.620,0.361]
select surf_pocket5, protein and id [2982,2983,2991,2996,3002,3006,3088,3261,5035,5020,3005,2999,5030,2976,2986,3138,3139,3140,5040,5043,5046,5047,5049,5050,3111,3106,3107,3108,3119,3120,3576,3244,3560,5288,5152,5155,5163,5151,5620,5164,5305,5604,5182,5183,5184] 
set surface_color,  pcol5, surf_pocket5 
   

deselect

orient
