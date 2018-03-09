

using Base.Test

@testset "HWintegration Unit Tests" begin
	@testset "testing components" begin
		@testset "check demand" begin
			@test 1==1
		end


		@testset "check gauss adjustment" begin
			@test HWintegration.coef1 == 3/2
			@test HWintegration.coef2 == 5/2
		end

		@testset "eqm condition for Q2" begin

		end
	end

	@testset "Testing Results" begin
	end

end
