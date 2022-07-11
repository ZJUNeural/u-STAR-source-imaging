function [Patch,pweight] = PatchGenerate(Seed,VertConn,VertArea,AreaDef)
numiter = 1;
pweight = [Seed,numiter];
Patch = Seed;
Area = sum(VertArea(Patch));
while Area <= AreaDef
    numiter = numiter+1;
    if Area <= AreaDef
        newverts = tess_scout_swell(Patch, VertConn);
        if ~isempty(newverts)
        Nouter = union(Patch,newverts);
        else
            break;
        end
        Area = sum(VertArea(Nouter));
    end
    if Area > AreaDef
        Ndiff = setdiff(Nouter,Patch);
        for i = 1: numel(Ndiff)
            Patch = union(Patch,Ndiff(i));
            nlist = [Ndiff(i),numiter];
            pweight = [pweight;nlist];
            Area = sum(VertArea(Patch));
            if Area > AreaDef
                break;
             end
        end
    else
         Patch = Nouter;
         nlist = [newverts(:),numiter*ones(length(newverts),1)];
         pweight = [pweight;nlist];
    end  
end
