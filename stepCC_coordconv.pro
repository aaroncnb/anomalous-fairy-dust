AME = READ_CSV('../Data/AME.txt', COUNT = amy, HEADER = amyhead, MISSING_VALUE='

glon = AME.field02[cerberus]
glat = AME.field03[cerberus]

AI = glon
BI = glat


EULER, AI, BI, AO, BO, [SELECT = 6 ]
