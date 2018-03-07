
module HWintegration

	const A_SOL = 4

	# imports
	using FastGaussQuadrature
	using Roots
	using Sobol
	using Plots
	using Distributions
	using NullableArrays
	using DataStructures  # OrderedDict


	# set random seed
	srand(12345)

	# demand function
	q(p) = 2p^(-1/2)

	# makes a plot for questions 1
	function plot_q1()
		plot(q,0.5,5,label="q(p)", title="Demand curve", xaxis="p", yaxis="quantity")
		hline!([1, 2], label=["qty at p= 4", "qty at p'= 1"], color = ["red", "green"])
	end

	function question_1b(n)
		#compute node points and weights
		x = gausslegendre(n)
		#Change function interval
		nodes_GL = 3/2 * x[1] + 5/2
		values_GL = [q(i) for i in nodes_GL]
		#Compute approximate integral
		GL_int = 3/2 * sum([x[2][i] * values_GL[i] for i in 1:n])
		deviation = A_SOL - GL_int
		# Plot graph
		global GL_graph = scatter(nodes_GL, values_GL, title="Gauss Legendre, n=$n")
		println("Gauss-Legendre method with $n nodes yields a change of $GL_int.
		Deviation from analytical solution is $deviation")
	end

	function question_1c(n)
		nodes_MC = rand(Uniform(1,4), n)
		values_MC = [q(x) for x in nodes_MC]
		MC_int = 3/n * sum(values_MC)
		deviation = A_SOL - MC_int

		#graph
		global MC_graph = scatter(nodes_MC, values_MC, title="Monte Carlo, n=$n")

		println("Monte Carlo method with $n nodes yields a change of $MC_int.
		Deviation from analytical solution is $deviation")
	end

	function question_1d(n)
		# Generate n points in one dimension (rescaled to be evenly spread on [1,4] interval)
		s = SobolSeq(1)
		nodes_qMC = 3*[[next(s) for i in 1:n][x][1] for x in 1:n] + 1
		values_qMC = [q(i) for i in nodes_qMC]
		qMC_int = 3/n * sum(values_qMC)
		deviation = A_SOL - qMC_int

		#graph
		global qMC_graph = scatter(nodes_qMC, values_qMC, title="Quasi Monte Carlo, n=$n")

		println("Quasi Monte Carlo method with $n nodes yields a change of $qMC_int.
		Deviation from analytical solution is $deviation")
	end

	# question 2

	# eqm condition for question 2
	# theta_1*p^(-1) + theta_2*p^(-1/2) = 2
	eqm(theta, supply, p) = theta[1]*p + theta[2]*p - supply

	function question_2a(n)
		mu = [0, 0]
		sigma = [0.02 0.01; 0.01 0.01]

		# Gauss-Hermite nodes
		#nodes = Any[]
		#push!(nodes,repeat(gausshermite(n),inner=[1],outer=[n]))  # dim1
		#push!(nodes,repeat(rules["hermite"][1],inner=[3],outer=[3]))  # dim2
		#weights = kron(rules["hermite"][2],kron(rules["hermite"][2],rules["hermite"][2]))
		#df = hcat(DataFrame(weights=weights),DataFrame(nodes,[:dim1,:dim2]))

		println("I could not solve this question")
	end

	function question_2b(n)
		println("I could not solve this question")
	end

	function question_2bonus(n)
		println("I could not solve this question")
	end

	# function to run all questions
	function runall()
		info("Running all of HWintegration")
		info("question 1:")
		plot_q1()
		for n in (10,15,20)
			info("============================")
			info("now showing results for n=$n")
			info("question 1b:")
			question_1b(n)	# make sure your function prints some kind of result!
			info("question 1c:")
			question_1c(n)
			info("question 1d:")
			question_1d(n)
			println("")
			info("question 2a:")
			q2 = question_2a(n)
			println(q2)
			info("question 2b:")
			q2b = question_2b(n)
			println(q2b)
			info("bonus question: Quasi monte carlo:")
			q2bo = question_2bonus(n)
			println(q2bo)
		end
	end
	info("end of HWintegration")

end
