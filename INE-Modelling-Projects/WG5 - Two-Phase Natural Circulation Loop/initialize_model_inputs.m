clear all
clc

p = [5: 2.5 :150];
Q = [13, 23, 33, 43];

inputs = combvec(p,Q);

for i=1:size(inputs,2)
    inputs(3,i) = find_FR_theory(inputs(1,i),inputs(2,i)*1000);
end

input_to_model = real(inputs);

save('inputs_to_model.mat', 'input_to_model');