
module HWintegration

	const A_SOL = 4  # enter your analytic solution.

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
	p=Any[]
	# demand function
    q(x)=2/(sqrt(x))
	a=3/2
	b=5/2
	g(x)=a*q(a*x+b)

    # Question 1

	# gauss-legendre adjustment factors
	a=3/2
	b=5/2
	# eqm condition for question 2

	# makes a plot for questions 1
	function plot_q1()
		plot(q,0,4,xlim=(0,4), ylim=(0,4),label="Demand")
		hline!([q(1)],label="p=1")
		hline!([q(4)],label="p=2")
	end


	function question_1b(n)
		println(sum(gausslegendre(n)[2][i]*g(gausslegendre(n)[1][i]) for i in 1:n )),
	    push!(p,plot(g,-1,1,xlim=(-1,1), ylim=(0,4),title="GaussLegendre",legend=false))
	    z=[]
	    for i in 1:n
	        push!(z,g(gausslegendre(n)[1][i]))
	    end
	    scatter!(gausslegendre(n)[1],z,legend=false)
	end


	function question_1c(n)
		MonteCarlo=rand(Uniform(1,4),n)
	    println(3*sum((1/n)*q(MonteCarlo[i]) for i in 1:n )),
	    push!(p,plot(q,0,4,xlim=(0,4), ylim=(0,4),title="MonteCarlo",legend=false))
	    z=[]
	    for i in 1:n
	        push!(z,q(MonteCarlo[i]))
	    end
	    scatter!(MonteCarlo,z,legend=false)

	end

	function question_1d(n)
        QMC=ScaledSobolSeq(1,[1.0],[4.0])
        QuasiMonteCarlo=([next(QMC) for i in 1:n])
        QMCX=[QuasiMonteCarlo[i][1] for i in 1:n]
        println(3*sum((1/n)*q(QMCX[i]) for i in 1:n )),
        push!(p,plot(q,0,4,xlim=(0,4), ylim=(0,4),title="QuasiMonteCarlo",legend=false))
        z=[]
        for i in 1:n
            push!(z,q(QuasiMonteCarlo[i][1]))
        end
        scatter!(QMCX,z,legend=false)


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
		plot_q1
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
			info("Creating Plot")
		end
		plot(p...)
	end
	info("end of HWintegration")

end
