empty_cells = c(
1:7, 10:27,
28:34, 40, 43:54,
55:61, 71:81,
82:88, 99:108,
109:115, 127:135,
136:142, 155:162,
163:168, 182:189,
190:195, 209:216,
217:221, 236:243,
244:248, 264:270,
271:273, 292:297, 
298:299, 320:324,
325:326, 350:351,
352, 377,
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
	378, #NS3
	375, #NB2
	359, #SD2
	370, #QB3
	343, #ON1
	344, #VT1
	348, #NS2
	280, #IA2
	276, #NE1
	255, #IL1
	258, #OH1
	227, #MO1
	234, #VA2
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
	414, #MB1
	284) #Michigan (Russ forest)

sampn = c(
	26, #MB2
	26, #AB1
	24, #ND2
	30, #NS3
	25, #NB2
	29, #SD2
	30, #QB3
	28, #ON1
	26, #VT1
	24, #NS2
	30, #IA2
	30, #NE1
	28, #IL1
	32, #OH1
	29, #MO1
	27, #VA2
	28, #KS1
	29, #VA1
	29, #TN1
	27, #SC1
	27, #MS2
	29, #LA1
	30, #TX2
	30, #QB1
	27, #WI2
	27, #WY2
	22, #MB1
	4)  #Michigan (Russ forest)

pop_name = c(	
	"MB2",
	"AB1",
	"ND2",
	"NS3",
	"NB2",
	"SD2",
	"QB3",
	"ON1",
	"VT1",
	"NS2",
	"IA2",
	"NE1",
	"IL1",
	"OH1",
	"MO1",
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
	"MB1",
	"Michigan")

popDF = data.frame(id = pop_name, grid.cell = samp_pops, sample.size = sampn)
popDF = popDF[order(popDF$grid.cell, decreasing = FALSE),]
