empty_cells = c(
1:7, 10:27,
28:34, 40, 43:54,
55:61, 71:81,
82:88, 99:108,
109:115, 127:135,
136:142, 155:162,
163:168, 182:189,
190:195, 209:216,
217:221,
	227, #MO1
236:243,
244:248, 264:270,
271:273,
284, #Michigan
292:297, 
298:299, 320:324,
325:326,
	348, #NS2
350:351,
352, 377,
	378, #NS3
379, 390:398, 403:405,
418:425, 427, 431:432,
442:459,
469:486,
495:513,
514:517, 521:540)

samp_pops = c(
	440, #MB2
	433, #AB1
	384, #ND2
#	378, #NS3
	375, #NB2
	359, #SD2
	370, #QB3
	343, #ON1
	344, #VT1
#	348, #NS2
	280, #IA2
	276, #NE1
	255, #IL1
	258, #OH1
	234, #VA2
#	227, #MO1
	225, #KS1
	206, #VA1
	150, #TN1
	98,  #SC1
	92,  #MS2
	38,  #LA1
	9,   #TX2
	426, #QB1
	363, #WI2
	328, #WY2
	414) #MB1
#	284) #Michigan (Russ forest)



    sampn = c(
	15, #MB2
	15, #AB1
#	15, #ND2
	15, #NS3
	15, #NB2
	14, #SD2
	15, #QB3
	15, #ON1
	15, #VT1
#	24, #NS2
	15, #IA2
	15, #NE1
	15, #IL1
	15, #OH1
#	1, #MO1
	15, #VA2
	15, #KS1
	15, #VA1
	15, #TN1
	15, #SC1
	15, #MS2
	15, #LA1
	15, #TX2
	15, #QB1
	15, #WI2
	15, #WY2
	15) #MB1
#	2   #Michigan (Russ forest)
 

pop_name = c(	
	"MB2",
	"AB1",
	"ND2",
#	"NS3",
	"NB2",
	"SD2",
	"QB3",
	"ON1",
	"VT1",
#	"NS2",
	"IA2",
	"NE1",
	"IL1",
	"OH1",
#	"MO1",
	"VA2",
	"KS1",
	"VA1",
	"TN1",
	"SC1",
	"MS2",
	"LA1",
	"TX2",
	"QB1",
	"WI2",
	"WY2",
	"MB1")
#	"Michigan")

popDF = data.frame(id = pop_name, grid.cell = samp_pops, sample.size = sampn)
popDF = popDF[order(popDF$grid.cell, decreasing = FALSE),]
