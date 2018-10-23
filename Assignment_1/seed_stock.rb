#The class Stockdatabase has 5 attributes and the function 
class StockDatabase
     @@number_of_stock = 0
    attr_accessor :seed_Stock  
    attr_accessor :mutant_Gene_ID
    attr_accessor :last_Planted
    attr_accessor :storange
    attr_accessor :grams_Remaining
    def initialize params={}
        @@number_of_stock+=1
        @seed_Stock = params.fetch(:seed_Stock, 'seed_Stock_unknow')
        @mutant_Gene_ID = params.fetch(:mutant_Gene_ID, 'mutant_Gene_id_unknow')
        @last_Planted = params.fetch(:last_Planted, 'unknow date')
        @storange = params.fetch(:storange, 'unknow storange')
        @grams_Remaining = params.fetch(:grams_Remaining, 'unknow remaining grams')
    end
    def how_many 
      return @@number_of_stock
    end
    
    #We create the function planting to simulate the loss in the stock of 7 grams of seeds plated
    def planting 
        grams_left=grams_Remaining.to_i-(7)
        if grams_left<=0
            return "#{seed_Stock}"
        else
            return "not true"
        end
    end
    
end
