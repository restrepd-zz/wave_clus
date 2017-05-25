function permuted_out=permut_part(to_permute)

num_to_display=500;
if to_permute>num_to_display
    permuted=randperm(to_permute);
    permuted_out=permuted(1:num_to_display);
else
    permuted_out=randperm(to_permute);
end
