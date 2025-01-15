clear all; clf

p = linspace(0, 1);
b = [p; 1-p];
output = dlmread('output.dat', ",", 0, 1);
output_data = output(:,1:2);

for i = 1:size(output_data,1)
    hold on; pp=plot(p,output_data(i,:)*b,'b');
end
box on
xlabel('Belief State')
ylabel('Cost-To-Go')



