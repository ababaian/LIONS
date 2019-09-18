#! /bin/awk -f

BEGIN{
  FS=OFS="\t";
  t = 0
}; 

{ 
  s = $1; 
  if (s in b){
    t++
  } else {
    t=1
  }; 
  b[s] = t;
  $1 = $1"_"b[s]; 
  a[$1] = $0; 
}; 

END{
  for (i in a) {
      print a[i]
  }
};
