
module HWintegration

	const A_SOL = -4

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

    q(x)=2/sqrt(x)

	# gauss-legendre adjustment factors for map change

    a = 4
    b = 1
    coef1 = (a-b)/2
    coef2 = (a+b)/2

	# eqm condition for question 2

	# makes a plot for questions 1
	function plot_q1()

	end


	function question_1b(n)

	end


	function question_1c(n)


	end

	function question_1d(n)



	end

	# question 2

	function question_2a(n)



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
			# info("question 1b:")
			# question_1b(n)	# make sure your function prints some kind of result!
			# info("question 1c:")
			# question_1c(n)
			# info("question 1d:")
			# question_1d(n)
			# println("")
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
