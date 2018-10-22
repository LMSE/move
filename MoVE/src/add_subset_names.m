function sys= add_subset_names(sys, react_name)

sys.sub_name= cell(size(sys.sub, 1), 1);
for i= 1:size(sys.sub, 1)
  sys.sub_name{i}= '{ ';
  for j= find(sys.sub(i, :))
    sys.sub_name{i}= [sys.sub_name{i}, sprintf('%g %s(%d) ', full(sys.sub(i, j)), react_name{j}, j)];
  end
  sys.sub_name{i}= [sys.sub_name{i}, '}'];
end
sys.react_name= react_name;
