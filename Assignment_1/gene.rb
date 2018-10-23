#Gene class
class Gene
    @@number_of_gene = 0
    attr_accessor :gene_ID 
    attr_accessor :gene_name
    attr_accessor :mutant_phenotype
    def initialize params={}
        @@number_of_gene+=1
        @gene_ID = params.fetch(:gene_ID, 'unknow gene ID')
        @gene_name = params.fetch(:gene_name, 'unknow gene name')
        @mutant_phenotype = params.fetch(:mutant_phenotype, 'unknow mutant phenotype')
    end
    def how_many 
      return @@number_of_gene
    end
    
    #def test_ID ()
      #if/A[Tt]\d[Gg]\d\d\d\d\d/ =~ ()
         # puts "Correct Gene_ID"
      #else
         # puts "Wrong Gene_ID"
      #end
    #end
end