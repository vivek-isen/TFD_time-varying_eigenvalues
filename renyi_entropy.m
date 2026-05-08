function rentropy = renyi_entropy(TFD,alpha)
TFD1 = TFD.^alpha;
num = sum(sum(TFD1));
den = (sum(sum(TFD)));
rentropy = log2(num/den)/(1-alpha);
end