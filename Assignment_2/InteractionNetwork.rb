#Here we define the Interaction Network Class
class InteractionNetwork
    @@numbers_of_genes = 0
    attr_accessor :interaction_IDs 
    def initialize params={}
        @@numbers_of_genes+=1
        @interaction_IDs = params.fetch(:interaction_IDs, 'unknown genes IDs')
    end
    def how_many 
      return @@numbers_of_genes
    end
end