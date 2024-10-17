import LibTiePie.Const.*
import LibTiePie.Enum.*

if ~exist('LibTiePie', 'var')
  % Open LibTiePie:
  LibTiePie = LibTiePie.Library;
else
    disp('Library Connection Failed');
end
disp('Done')

%%

test()

function test()
    import LibTiePie.Const.*
    import LibTiePie.Enum.*

    LibTiePie = LibTiePie.Library;
disp('Done')

end