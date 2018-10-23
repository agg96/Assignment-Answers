#The class Cross has 6 attributes and the chi-square function 
class Cross
    @@number_of_cross = 0
    attr_accessor :parent1 
    attr_accessor :parent2
    attr_accessor :f2_wild
    attr_accessor :f2_P1
    attr_accessor :f2_P2
    attr_accessor :f2_P1P2
    def initialize params={}
        @@number_of_cross+=1
        @parent1 = params.fetch(:parent1, 'unknown parent1')
        @parent2 = params.fetch(:parent2, 'unknown parent2')
        @f2_wild = params.fetch(:f2_wild, 'unknown f2_wild')
        @f2_P1 = params.fetch(:f2_P1, 'unknown f2_P1')
        @f2_P2= params.fetch(:f2_P2, 'unknown f2_P2')
        @f2_P1P2 = params.fetch(:f2_P1P2, 'unknown f2_P1P2')
    end
    def how_many 
      return @@number_of_cross
    end
    # We define the chi-square function. First we calculate the expected data and then we use the chi-square formula
    def chi_square
      total_seeds = f2_wild.to_f + f2_P1.to_f + f2_P2.to_f + f2_P1P2.to_f
      f2_wild_expect = total_seeds*9/16
      f2_P1_expect = total_seeds*3/16
      f2_P2_expect = total_seeds*3/16
      f2_P1P2_expect = total_seeds*1/16
      result = ((f2_wild.to_f - f2_wild_expect)**2 / f2_wild_expect) + ((f2_P1.to_f - f2_P1_expect)**2 / f2_P1_expect) + ((f2_P2.to_f - f2_P2_expect)**2 / f2_P2_expect) +((f2_P1P2.to_f - f2_P1P2_expect)**2 / f2_P1P2_expect)
      return result
    end
end