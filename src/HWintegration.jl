
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
        q_trans(x) = coef1 * 2/sqrt(coef1*x+coef2)
        rule = gausslegendre(n)
        nodes = values(rule)[1,]
        weights = values(rule)[2,]
        weighted_q = weights .* map(q_trans, nodes)
        approx = sum(weighted_q)
        println("The approximated change in consumer surplus when n=$n is $approx.")
        println("The distance between the true result and the approximation when n=$n is ", 4 - approx)
        scatter(nodes, map(q_trans, nodes), label = "Gauss-Legendre", xlab ="Integration nodes", ylab="Function Value" )
        plot!(q_trans, label = "Scaled Demand Function")
    end


    function question_1c(n)

        function rand_uniform(a, b, n)
            a + rand(n)*(b - a)
        end

        q_trans(x) = coef1 * 2/sqrt(coef1*x+coef2)
        random = rand_uniform(-1, 1, n)
        x = sort(random)
        y = map(q_trans, x)
        approx = sum(y)/n
        println("Results using Monte Carlo integration : ")
        println("The approximated change in consumer surplus when n=$n is $approx.")
        println("The distance between the true result and the approximation when n=$n is ", 4 - approx)
        scatter!(x, y, label = "Monte Carlo", xlab ="Integration nodes", ylab="Function Value" )

    end


    function question_1d(n)

        #= We define a function Sobol that generates the sobol
        sequence in a usable way =#

        function sobol(j)
           s = ScaledSobolSeq(1, [-1.0], [1.0])
           seq = hcat([next(s) for i in 1:j]...)
           a = []
            for i in 1:j
                push!(a, seq[1, i])
            end
        return a
        end

        q_trans(x) = coef1 * 2/sqrt(coef1*x+coef2)

        rand = sobol(n)
        x = sort(rand)
        y = map(q_trans, x)

        approx = sum(y)/n
        println("Results using Quasi Monte Carlo integration : ")
        println("The approximated change in consumer surplus when n=$n is $approx.")
        println("The distance between the true result and the approximation when n=$n is ", 4 - approx)
        scatter!(x, y, label = "Q. Monte Carlo", xlab ="Integration nodes", ylab="Function Value" )

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
        info("question 2:")
        plotlyjs()
        plots = []
        for n in (10,15,20)
            info("============================")
            info("now showing results for n=$n")
            info("question 1b:")
            question_1b(n)
            info("question 1c:")
            question_1c(n)
            info("question 1d:")
            push!(plots, question_1d(n))
        end
        plot(plots[1], plots[2], plots[3])


            info("question 2:")
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
