
module HWintegration

	const A_SOL = 4  # enter your analytic solution. -18 is wrong.

	# imports
	using FastGaussQuadrature
	using Roots
	using Sobol
	using Plots
	using NullableArrays
	using DataStructures  # OrderedDict
	using Distributions


	# set random seed
	srand(12345)

	# price shocks q1
	a = 1
	b = 4
	#supply for q2
	supply = 2
	#variance covariance matrix for q2
	sigma = [[0.02, 0.01] [0.01, 0.01]]
	# demand function
	q(p) = 2./sqrt.(p)
	# gauss-legendre adjustment factors for map change
	b_a_2 = (b-a)/2
    ab_2 = (a+b)/2
	# eqm condition for question 2
	eqm(theta, supply, p) = theta[1]*p + theta[2]*p - supply
	# makes a plot for questions 1
	function plot_q1()
		p_grid = vcat([p  for p in 0.1:0.1:6]...)
		demand = q(p_grid)
		plot(p_grid, demand, xaxis="Price", yaxis="Demand", label="q(d)", lw = 3)
		plot!(4*ones(Int32, size(p_grid,1)), p_grid, label="equilibrium price 4", lw=2, ls=:dash, color=:orange )
		plot!(ones(Int32, size(p_grid,1)), p_grid, label="equilibrium price 1", lw=2, ls=:dash, color=:orange )
	end


	function question_1b(np) #Quadrature by Gauss Legendre
		w_sum_gl = Array{Float64}(size(np,2),1) #initializing array to store surplus values
		for i in 1:size(np,2)
			glw = collect(gausslegendre(np[i])) #transform tuple structure into array structure
			glw_nodes = glw[1] #isolate nodes array
			glw_weights = glw[2] #isolate weights array
			demand_app = q.((b_a_2*glw_nodes+ab_2)) #broadcast demand function over properly scaled nodes, transform is needed since GL works only for integrals between -1 and 1
			w_sum_gl[i] = b_a_2*glw_weights.'*demand_app #compute weighted sum by matrix multiplication. Scale factor is result of the interval transformation in order to use GL method
			info("for $(np[i]) number of points the approx is:")
			info(w_sum_gl[i])
		end
		plot(w_sum_gl)
	end


	function question_1c(np) #Montecarlo
		w_sum_mc = Array{Float64}(size(np,2),1) #initializing array to store surplus values
	    for i in 1:size(np,2)
	        rand_p_01 = rand(np[i]) #np random numbers
	        #We use a uniform distribution to draw points in the interval (a,b). Recall uniform cdf = (x - a)/(b - a)
	        rand_p_ab = (b-a)*rand_p_01 + a #Applying a quantile formula for the uniform distribution i our sample space (a,b)
	        w_sum_mc[i] = (b-a)*sum(q.(rand_p_ab))/np[i] #weighting the sum for 1/N and the volume of the sample space (b-a)
			info("for $(np[i]) number of points the approx is:")
			info(w_sum_mc[i])
	    end
		plot(w_sum_mc)
	end

	function question_1d(np)
		w_sum_qmc = Array{Float64}(size(np,2),1) #initializing array to store surplus values
	    sb = SobolSeq(1) #Obtaining a 1 dimensional sobol sequence
	    for i in 1:size(np,2)
	        rand_p_sobol = 0
	        rand_p_sobol = vcat([next(sb) for i = 1:np[i]]...) #building an array of np elements from the sobol sequence
	        rand_p_ab = (b-a)*rand_p_sobol + a #Applying a quantile formula taking elements of sobol seq into a uniform distribution i our sample space (a,b)
	        w_sum_qmc[i] = (b-a)*sum(q.(rand_p_ab))/np[i] #weighting the sum for 1/N and the volume of the sample space (b-a)
			info("for $(np[i]) number of points the approx is:")
			info(w_sum_qmc[i])
	    end
		plot(w_sum_qmc)
	end

	# question 2

	function question_2a(n)
		unt_theta1 = repeat(gausshermite(n)[1], inner=[1], outer=[n]) 	# repeating gausshermite to get 1 dimension of nodes
		unt_theta2 = repeat(gausshermite(n)[1], inner=[n], outer=[1]) 	# repeating gausshermite to get 2nd dimension of nodes
		unt_theta = [unt_theta1 unt_theta2] 	#concatanating both dimensions of nodes into single array
		trf_log_theta = chol(sigma)*unt_theta.' 	# logaritmic versions of transformed thetas
		trf_theta = exp.(trf_log_theta).' 	#true value of thetas
		gh_weights = kron(gausshermite(n)[2],gausshermite(n)[2]) 	#obtaining weights form gausshermite and applying directly kronecker product
		p_nodes = Array{Float64}(n^2,1)	#initializing prices nodes for later computation
		for i in 1:n^2
		    p_nodes[i] = fzero(p -> eqm(trf_theta[i,:], supply, p),1) #finding value of price for each pair of realizations of thetas with given supply
		end
		exp_price_gh = gh_weights.'*p_nodes	 #computation of expectation
		exp_var_gh = gh_weights.'*((p_nodes - mean(p_nodes)).^2) 	# computation of variance
		info("the expectation of the prices is $exp_price_gh")
		info("the variance of the prices is $exp_var_gh")
	end

	function question_2b(n)
		unt_theta = randn(n^2,2)
		trf_log_theta = chol(sigma)*unt_theta.' 	# logaritmic versions of transformed thetas
		trf_theta = exp.(trf_log_theta).'
		p_nodes = Array{Float64}(n^2,1)	#initializing prices nodes for later computation
		for i in 1:n^2
		    p_nodes[i] = fzero(p -> eqm(trf_theta[i,:], supply, p),1) #finding value of price for each pair of realizations of thetas with given supply
		end
		exp_price_mc = 1/(n^2)*sum(p_nodes) #computation of expectation
		exp_var_mc = 1/(n^2)*sum((p_nodes - mean(p_nodes)).^2) # computation of variance
		info("the expectation of the prices is $exp_price_mc")
		info("the variance of the prices is $exp_var_mc")
	end

	function question_2bonus(n)

	end

	# function to run all questions
	function runall()
		# Number of point to test functions
		np = [10 15 20]
		# Showing results
		info("Running all of HWintegration")
		info("question 1:")
		plot_q1()
		info("============================")
		info("Integrating by Quadrature using GL rule")
		question_1b(np)
		info("============================")
		info("Integrating by Monte Carlo")
		question_1c(np)
		info("============================")
		info("Integrating by Quasi Monte Carlo")
		question_1d(np)
		for n in (10,15,20)
			info("============================")
			info("now showing results for n=$n")
			info("question 2a:")
			question_2a(n)
			info("question 2b:")
			question_2b(n)
			info("bonus question: Quasi monte carlo:")
			question_2bonus(n)
			println("I could not solve this question")
		end
	end
	info("end of HWintegration")

end
