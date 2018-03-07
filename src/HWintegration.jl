
module HWintegration

	const A_SOL = 4

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

    q(x)=2/sqrt(x)

	# gauss-legendre adjustment factors for map change

    a = 4
    b = 1
    coef1 = (a-b)/2
    coef2 = (a+b)/2

	# eqm condition for question 2

    # S = X + D
    # theta_1/x + theta_2/sqrt(x) - 2 = 0

	# makes a plot for questions 1
	function plot_q1()
        gr()
        x=0.5:10
        plot(x, 2 ./ sqrt.(x))
        vline!([4], color=:green,lw=1,label="p = 4")
        vline!([1], color=:orange,lw=1,label="p = 1")
	end


    function question_1b(n)
        q_trans(x) = coef1 * 2/sqrt(coef1*x+coef2)
        rule = gausslegendre(n)
        nodes = values(rule)[1,]
        weights = values(rule)[2,]
        weighted_q = weights .* map(q_trans, nodes)
        approx = sum(weighted_q)
        println("When n = $n :")
        println("The approximated change in consumer surplus is $approx.")
        println("The distance between the true result and the approximation is ", A_SOL - approx)
        scatter(nodes, map(q_trans, nodes), label = "Gauss-Legendre", xlab ="Integration nodes", ylab="Function Value" )
        plot!(q_trans, label = "Scaled Demand Function")
    end


    function question_1c(n)

        q_trans(x) = coef1 * 2/sqrt(coef1*x+coef2)
        x = sort(rand(n))
        y = map(q_trans, x)
        approx = sum(y)*(1+1)/n
        println("When n = $n :")
        println("The approximated change in consumer surplus is ", approx)
        println("The distance between the true result and the approximation is ", A_SOL - approx)
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

        approx = sum(y)*(1+1)/n
        println("When n = $n :")
        println("The approximated change in consumer surplus is $approx.")
        println("The distance between the true result and the approximation is ", A_SOL - approx)
        scatter!(x, y, label = "Q. Monte Carlo", xlab ="Integration nodes", ylab="Function Value" )

    end

	# question 2

	function question_2a(n)
        rule = gausshermite(n)
        nodes = Any[]
        push!(nodes,repeat(rule[1],inner=[1],outer=[n]))
        push!(nodes,repeat(rule[1],inner=[n],outer=[1]))
        #Making the nodes being of variance 1
        nodes[1] = nodes[1]/sqrt(var(nodes[1]))
        nodes[2] = nodes[2]/sqrt(var(nodes[2]))
        weights = kron(rule[2],rule[2])
        points = DataFrame(weights = weights)
        points[:dim1] = nodes[1]
        points[:dim2] = nodes[2]
        theta_unnorm = hcat(points[2], points[3])

        #= Making log(theta1) and log(theta2) such that
        var(log(theta1)) = 0.01
        var(log(theta2)) = 0.02
        cov(log(theta1), log(theta2)) = 0.01 =#
        dim = n^2
        SIG = [sqrt(0.01) 0 ; sqrt(0.01) sqrt(0.01)]
        theta = (SIG*theta_unnorm')'
        points[:theta_1] = exp.(theta[1:dim,])
        points[:theta_2] = exp.(theta[dim+1:dim*2,])

        # Generating price vector

        price = []
        for i in 1:dim
            a = points[:theta_1][i]
            b = points[:theta_2][i]
            push!(price, fzero(x -> a/x + b/sqrt(x) - 2, 0, 10))
        end

        points[:price] = price

        approx_E = sum(points[:price] .* points[:weights])
        approx_V = sum((points[:price] .^ 2) .* points[:weights])
        println("When n = $n :")
        println("The approximated expected price is $approx_E.")
        println("The approximated expected variance is $approx_V")
	end

	function question_2b(n)
        L = 0.1 * [1/sqrt(2) 0; 1/sqrt(2) sqrt(2)]
        theta_1 = randn(n)
        theta_2 = randn(n)
        theta_unnorm = hcat(theta_1, theta_2)
        theta_norm = (L*theta_unnorm')'
        theta_1_norm = [exp.(theta_norm[1:n,])]
        theta_2_norm = [exp.(theta_norm[n+1:n*2,])]
        points = DataFrame(theta_1_norm = theta_1_norm[1], theta_2_norm = theta_2_norm[1])

        price = []
        for i in 1:n
            a = points[:theta_1_norm][i]
            b = points[:theta_2_norm][i]
            if (a<0) && (b<0)
                push!(price,0)
            else
                push!(price, fzero(x -> a/x + b/sqrt(x) - 2, 0, 10))
            end
        end
        points[:price] = price
        approx_E = sum(points[:price])/n
        approx_V = sum( (points[:price]- approx_E) .^ 2 ) / n
        println("When n = $n :")
        println("The approximated expected price is $approx_E.")
        println("The approximated expected variance is $approx_V")

	end

	function question_2bonus(n)

        function sobol(j)
           s = ScaledSobolSeq(1, [-1.0], [1.0])
           seq = hcat([next(s) for i in 1:j]...)
           a = []
            for i in 1:j
                push!(a, seq[1, i])
            end
        return a
        end

        nodes = sobol(n)
        theta_unnorm = []
        push!(theta_unnorm,repeat(nodes,inner=[1],outer=[n]))
        push!(theta_unnorm,repeat(nodes,inner=[n],outer=[1]))
        L = 0.1 * [1/sqrt(2) 0; 1/sqrt(2) sqrt(2)]
        points = DataFrame(theta_1 = theta_unnorm[1], theta_2 = theta_unnorm[2])
        theta = hcat(points[1], points[2])
        theta_norm = (L*theta')'
        dim = n^2
        theta_1_norm = [theta_norm[1:dim,]]
        theta_2_norm = [theta_norm[dim+1:dim*2,]]
        points[:theta_1_norm] = exp.(theta_1_norm[1])
        points[:theta_2_norm] = exp.(theta_2_norm[1])
        price = []
            for i in 1:dim
                a = points[:theta_1_norm][i]
                b = points[:theta_2_norm][i]
                    if (a<0) && (b<0)
                        push!(price,0)
                    else
                        push!(price, fzero(x -> a/x + b/sqrt(x) - 2, 0, 10))
                    end
            end
        points[:price] = price
        approx_E = sum(points[:price])/dim
        approx_V = sum( (points[:price]- approx_E) .^ 2 ) / dim
        println("When n = $n :")
        println("The approximated expected price is $approx_E.")
        println("The approximated expected variance is $approx_V")
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
            info("question 1b (Gauss Legendre):")
            question_1b(n)
            info("question 1c (Monte Carlo):")
            question_1c(n)
            info("question 1d (Quasi Monte Carlo):")
            push!(plots, question_1d(n))

            info("question 2:")
			info("question 2a (Gauss Hermite)")
			q2 = question_2a(n)
			println(q2)
			info("question 2b (Monte Carlo):")
			q2b = question_2b(n)
			println(q2b)
			info("bonus question (Quasi monte carlo):")
			q2bo = question_2bonus(n)
			println(q2bo)
		end
    plot(plots[1], plots[2], plots[3])
	end
	info("end of HWintegration")

end
