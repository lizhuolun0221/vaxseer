import importlib
import json
import sys
from pathlib import Path
import subprocess
root=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(root/'vaxseer'))
modules=['torch','torchvision','torchaudio','pytorch_lightning','lightning_transformers','transformers','torchmetrics','Bio','numpy','scipy','sklearn','pandas','seaborn','jupyterlab','openpyxl','dropbox','models','data.data_modules','config']
result={'python':sys.version,'executable':sys.executable,'imports':{}}
for name in modules:
    try:
        module=importlib.import_module(name)
        result['imports'][name]={'ok':True,'version':getattr(module,'__version__',None)}
    except Exception as exc:
        result['imports'][name]={'ok':False,'error':repr(exc)}
import torch
x=torch.randn(4,3,requires_grad=True)
(x.square().sum()).backward()
result['cpu_autograd']=bool(torch.isfinite(x.grad).all())
result['cuda_available']=torch.cuda.is_available()
result['torch_cuda_build']=torch.version.cuda
p=subprocess.run([sys.executable,'-m','pip','check'],capture_output=True,text=True)
result['pip_check']={'ok':p.returncode==0,'output':p.stdout+p.stderr}
p=subprocess.run([str(Path(sys.executable).parent/'mmseqs'),'version'],capture_output=True,text=True)
result['mmseqs']={'ok':p.returncode==0,'version':p.stdout.strip(),'stderr':p.stderr}
print(json.dumps(result,indent=2))
sys.exit(not(all(x['ok'] for x in result['imports'].values()) and result['cpu_autograd'] and result['pip_check']['ok'] and result['mmseqs']['ok']))
