classdef materialProperties < matlab.mixin.Copyable
properties 
    rho
    mu 
    beta
    sigma
end
methods
    
    
    
     function Add(obj, name, value)
         switch(name)
             case 'rho'
                 obj.rho=value;
             case 'mu'
                 obj.mu=value;
             case 'sigma'
                 obj.sigma=value;
             case 'beta'
                 obj.beta=value;
             otherwise
                 
        
     end
     
    
    
    end
end