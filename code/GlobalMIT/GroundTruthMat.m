function truth = GroundTruthMat(csvfile, m, n)
truth = zeros(m,n);
vec = csvread(csvfile);
for point=vec'
    truth(point(1), point(2)) = point(3);
end