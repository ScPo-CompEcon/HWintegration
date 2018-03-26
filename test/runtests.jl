

using Base.Test

@testset "HWintegration Unit Tests" begin
	@testset "testing components" begin
		@testset "check demand" begin
			@test 1==1
		end


		@testset "check gauss adjustment" begin
			@test HWintegration.coef1 == 1.5
			@test HWintegration.coef2 == 2.5
		end

		@testset "eqm condition for Q2" begin
		@test 1==1
		end
	end

	@testset "Testing Results" begin
	@test 1 ==1
	end

end
