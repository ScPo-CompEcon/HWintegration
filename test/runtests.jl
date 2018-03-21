

using Base.Test

@testset "HWintegration Unit Tests" begin
	@testset "testing components" begin
		@testset "check demand" begin
			@test HWintegration.q(4) == 1
		end


		@testset "check gauss adjustment" begin
		end

		@testset "eqm condition for Q2" begin
			@test egm_cond.q(0,0,1) == 0
		end
	end

	@testset "Testing Results" begin
		@test HWintegration.question_1b(20) > 3.5
		@test HWintegration.question_1c(20) > 3.5
		@test HWintegration.question_1d(20) > 3.5
		@test HWintegration.question_1b(20) < 4.5
		@test HWintegration.question_1c(20) < 4.5
		@test HWintegration.question_1d(20) < 4.5
	end
    
end
