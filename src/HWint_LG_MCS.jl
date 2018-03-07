
module HWintegration

	const A_SOL = 4
	println("The analytical solution is $A_SOL")
	# imports
	using FastGaussQuadrature
	using Roots
	using Sobol
	using Plots
	using Distributions
	using NullableArrays
	using DataFrames
	using DataStructures  # OrderedDict

	# set random seed
	srand(12345)

	# demand function
	function q(p)
		2*p ^(-0.5)
	end

	# makes a plot for question 1
	function plot_q1()
	global dem = plot(q, 0.1, 5, title="Demand Function", xaxis=("Price"), yaxis=("Quantity"), labels = "q(p)")
	hline!([2, 1], label = "", linestyle=:dash)
	scatter!([1], [2], label = "", color = ["red"])
	scatter!([4], [1], label = "", color = ["red"])
	end

	function question_1a(n)
	nw = gausslegendre(n)
	glx = 3/2*nw[1]+ 5/2
	gly = [q(3/2*x + 5/2) for x in nw[1]]
	int = 3/2*sum([q(3/2*x + 5/2) for x in nw[1]] .* nw[2])
	global gl = scatter(glx, gly, labels = "Gauss-Legendre", title = "n = $n")
	diff = A_SOL-int
	println("The Gauss-Legendre method gives $int when using $n nodes. The error is $diff")
	end

	function question_1b(n)
	mcx = rand(Uniform(1,4), n)
	mcy = q.(mcx)
	int = 3/n * sum([q(x) for x in mcx])
	global mc = scatter(mcx, mcy, labels = "Monte-Carlo", title = "n = $n")
	println("The Monte-Carlo method gives $int")
	end

	function question_1c(n)
	s = SobolSeq(1)
	qmcx = 3*[ hcat([next(s) for i = 1:n])[x][1] for x = 1:n]+1
	qmcy = q.(qmcx)
	int = 3/n * sum([q(x) for x in qmcx])
	global qmc = scatter(qmcx, qmcy, labels = "Quasi Monte-Carlo", title = "n = $n")
	println("The quasi-Monte-Carlo method gives $int")
	end

	# question 2

	function question_2a()
		function p(t = [0,0])
			function eqm(x)
				exp(t[1])* float(x)^(-1) + exp(t[2])*float(x)^(-0.5) - 2
			end
			fzero(eqm, [0,10000])
		end
		μ = [0,0]
		Σ = [0.02   0.01 ; 0.01  0.02]
		Ω = chol(Σ)
		#rand(MvNormal(μ, Σ), n)
		nodes = Any[]
		push!(nodes,repeat(gausshermite(10)[1],inner=[1],outer=[10]))  # dim1
		push!(nodes,repeat(gausshermite(10)[1],inner=[2],outer=[5]))  # dim2
		#weights1 = kron(gausshermite(10)[2], gausshermite(10)[2])
		dict = Dict(gausshermite(10)[1][i] => gausshermite(10)[2][i] for i in 1:10)
		#df = hcat(DataFrame(weights=weights),DataFrame(nodes,[:dim1,:dim2]))
		coor = [ [nodes[1][i], nodes[2][i]] for  i in 1: 100]
		weights2 = [dict[coor[i][1]]*dict[coor[i][2]] for i in 1:length(coor)]
		cofv = [Ω *sqrt(2)*coor[i]+ μ for i in 1:100]
		expectation = (1/sqrt(π)  * sum(weights2 .* p.( cofv)))
		variance = (1/sqrt(π)  * sum(weights2 .* ((p.(cofv)- expectation).^2)))
		println("The Gauss-Hermite method for approximating integrals gives:")
		println("For the expectation: $expectation")
		println("For the variance: $variance")
	end

	function question_2b(n)
		μ = [0,0]
		Σ = [0.02   0.01 ; 0.01  0.02]
		randompoints = rand(MvNormal(μ, Σ),  n)
		randompoints = [[randompoints[1,i], randompoints[2,i]] for i in 1:n]
		function p(t = [0,0])
			function eqm(x)
				exp(t[1])* float(x)^(-1) + exp(t[2])*float(x)^(-0.5) - 2
			end
			fzero(eqm, [0,10000])
		end
		expectation = 1/n * sum(p.( randompoints))
		variance = 1/n * sum((p.(randompoints) - expectation).^2)
		println("The Monte-Carlo method for approximating integrals gives:")
		println("For the expectation: $expectation")
		println("For the variance: $variance")
	end

	function question_2bonus()
		println("We couldn't solve this question.")
	end

	# function to run all questions
	function runall()
		info("Running all of HWintegration")
		info("Question 1:")
		for n in (10,15,20)
			info("Now showing results for n = $n")
			info("Question 1.A - Gauss-Legendre:")
			question_1a(n)
			info("Question 1.B - Monte-Carlo")
			question_1b(n)
			info("Question 1.C - Quasi-Monte-Carlo")
			question_1c(n)
			println("")
			plot_q1()
			allplots = plot(dem, gl, mc, qmc, layout = (2,2))
			display(allplots)
			savefig("n=$n.png")
			println("Plot for n = $n added to the current directory as a .png file.")
		end
		info("============================")
		println("")
		info("Question 2.A - Gauss-Hermite:")
		question_2a()
		info("Question 2.B - Monte-Carlo:")
			for n in (10,15,20)
			info("Now showing results for n = $n")
			question_2b(n)
		end
		info("Bonus Question: Quasi Monte-Carlo:")
		question_2bonus()
	info("End of HWintegration")
	end
end
