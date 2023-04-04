	import numpy as np
	import therm

	#initializing models with guessed values for parameters
	th_h  = therm.EXI(1.0e10)
	th_q  = therm.pQCD(10.0, 5.0)
	th_SF = therm.Standard(150.0, 2)

	th = therm.Crossover(th_h, th_q, th_SF)

	rs = [2,3,4,5,6] #list of r values to test
	opitimal_list = [] #list of results for each r

	#get optimal parameters for each value of r
	for r in rs:
		th_SF.r_exponent=r
		optimal_list.append(therm.Optimize(th))

	#print results
	for i, r in enumerate(rs):
		print(r, "| chi^2:", optimal_list[i][0]," | ", optimal_list[i][1:])

	i_best = np.argmin(optimal_list[:][0])
	print("\nBest fit parameters")
	print(rs[i_best], optimal_list[i_best][1:])