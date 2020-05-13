module jWords

export words, slurpwords

function slurpwords(fin::String)
#Read a text file and get a count of distinct words.
  fid = open(fin)
  s = read(fid, String)  #slurps the whole file
  println(s)
  s = lowercase(s)
  close(fid)
  d = words(s)  #a dictionary of words with their counts
end

function words(s::String)
#Get the counts of distinct words in a string.
  index = firstindex(s)
  d = Dict()     #start a new dictionary with word-count entries

@label start
  if isletter(s[index])
    w = ""  #starts a new word
    @goto build_word
    #@goto new_word
  else
    if index == lastindex(s) @goto done end  #no more characters to process
    index = nextind(s, index)
    @goto start
  end

@label build_word
  w = w * s[index]
  if index == lastindex(s) @goto end_word end  #no more characters to process
  index = nextind(s, index)
  if isletter(s[index])
    @goto build_word
  else
    @goto end_word
  end

@label end_word
  if w ∈ keys(d)
    d[w] += 1  #increments the count value for this word (note square brackets)
  else
    d[w] = 1  #this word seen for the first time
  end
  if index == lastindex(s) @goto done end  #no more characters to process
  @goto start

@label done
return d

end  #endfunction
end  #endmodule
