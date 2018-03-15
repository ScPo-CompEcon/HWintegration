
module HWintegration

	const A_SOL = 4  # enter your analytic solution. -18 is wrong.

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

	q(p) = 2p^(-0.5)

	# gauss-legendre adjustment factors for map change

	a = 1
	b = 4
	adj1 = (b - a)/2
	adj2 = (a + b)/2

	# eqm condition for question 2

	# makes a plot for questions 1
	function plot_q1()
		plot(p->q(p),0.5,5,label="Q")
		vline!([1,4],color=:red)
	end


	function question_1b(n)
		nodes, weights = gausslegendre(n)
		nodes_points = Vector{Float64}(n)
		for i in 1:n
			nodes_points[i] = adj1*nodes[i] + adj2
			nodes[i] = q(adj1*nodes[i] + adj2)
		end
		plot_q1()
		vline!(nodes_points,color=:blue)
		res = adj1*dot(nodes,weights)
		println(res)
	end


	function question_1c(n)
		nodes_points = rand(n)*(b-a) + a
		nodes = Vector{Float64}(n)
		for i in 1:n
			nodes[i] = q(nodes_points[i])
		end
		plot_q1()
		vline!(nodes_points,color=:blue)
		res = dot(ones(n),nodes)*(b-a)/n
		println(res)
	end

	function question_1d(n)
		s = SobolSeq(1)
		nodes_points = hcat([next(s)*(b-a) + a for i = 1:n]...)'
		nodes = Vector{Float64}(n)
		for i in 1:n
			nodes[i] = q(nodes_points[i])
		end
		plot_q1()
		vline!(nodes_points,color=:blue)
		res = dot(ones(n),nodes)*(b-a)/n
		println(res)
	end

	# question 2

	function question_2a(n)
		repeat(gausshermite(n)[1], inner=[1], outer=[n]) 
		

	end

	function question_2b(n)

		
	end

	function question_2bonus(n)

		


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
			#info("question 2a:")
			#q2 = question_2a(n)
			#println(q2)
			#info("question 2b:")
			#q2b = question_2b(n)
			#println(q2b)
			#info("bonus question: Quasi monte carlo:")
			#q2bo = question_2bonus(n)
			#println(q2bo)
		end
	end
	info("end of HWintegration")

end
