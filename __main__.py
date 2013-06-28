#!/usr/bin/python
#coding=utf-8

from engine.simplex import LinearProblem, ProblemType, Operator, simplex

products = [["eggs", 157, 13, 12, 0],
            ["rice", 344, 6, 1, 79],
            ["apple", 47, 0.4, 0.4, 9.8],
            ["curd", 100, 18, 1.8, 3.3],
            ["milk", 35, 2.8, 0.5, 4.6],
            ["chicken", 100, 18, 1.8, 3.3],
            ["tuna", 101, 24, 1, 0],
            ["whey", 406, 75, 3, 15.6],
            ["buckwheat", 313, 12.6, 3.3, 62.1],
            ["cheese", 350, 23.5, 26, 0],]

def getFunction():
	func = []
	for i in xrange(len(products)):
		func.append(float(100))
	return func[:]

def getTerms():
    terms, calList, protList, fatList, carbList = [], [], [], [], []
    for i in xrange(len(products)):
        calList.append(float(products[i][1]))
        protList.append(float(products[i][2]))
        fatList.append(float(products[i][3]))
        carbList.append(float(products[i][4]))
    terms = [calList[:], calList[:], protList[:], protList[:], fatList[:], fatList[:], carbList[:], carbList[:]]
    return terms[:]
    
def getLimits(calFrom, calTo, protFrom, protTo, fatFrom, fatTo, carbFrom, carbTo):
	limits = [float(calFrom), float(calTo), float(protFrom), float(protTo), 
			  float(fatFrom), float(fatTo), float(carbFrom), float(calTo)]
	return limits[:]

def getOperators():
	operators = []
	for i in xrange(4):
		operators.append(Operator.enum['>='])
		operators.append(Operator.enum['<='])
	return operators[:]

def toResults(ans):
	result = []
	for i in xrange(len(products)):
		for j in xrange(len(ans[3])):
			if i == ans[3][j]:
				result.append([i, ans[2][j]])
	return result[:]

def printResults(result):
	calTotal, protTotal, fatTotal, carbTotal = 0, 0, 0, 0
	print "Products:"
	for i in xrange(len(result)):
		print products[result[i][0]][0], int(round(result[i][1]*100)), "grams"
		calTotal  += float(products[result[i][0]][1]) * result[i][1]
		protTotal += float(products[result[i][0]][2]) * result[i][1]
		fatTotal  += float(products[result[i][0]][3]) * result[i][1]
		carbTotal += float(products[result[i][0]][4]) * result[i][1]
	print "\nTotal:"
	print "calories", int(round(calTotal)), "kc"
	print "proteins", int(round(protTotal)), "grams"
	print "fats",     int(round(fatTotal)), "grams"
	print "carbs",    int(round(carbTotal)), "grams"

def addVariableTerm(productName, weigth, terms, lims, op, mode):
	op.append(Operator.enum[mode])
	lims.append(float(weigth))
	term = []
	for i in xrange(len(products)):
		if products[i][0] == productName:
			term.append(float(1))
		else:
			term.append(float(0));
	terms.append(term[:])

if __name__ == '__main__':
	terms = getTerms()
	func  = getFunction()
	lims  = getLimits(1900, 2100, 180, 200, 50, 60, 180, 200)
	op    = getOperators()

	addVariableTerm("milk", 10, terms, lims, op, '<=')
	addVariableTerm("milk", 5, terms, lims, op, '>=')

	addVariableTerm("eggs", 1.8, terms, lims, op, '=')

	addVariableTerm("whey", 0.64, terms, lims, op, '=')
	addVariableTerm("curd", 2, terms, lims, op, '=')
	
	lp = LinearProblem(ProblemType.enum['min'], func, terms, lims, op)
	res = simplex(lp)
	if res[0]:
		total = toResults(res)
		printResults(total)