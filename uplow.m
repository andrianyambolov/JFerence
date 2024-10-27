function [out] = uplow(instr)

out = char(instr);
out = string([upper(out(1)), out(2:end)]);

end