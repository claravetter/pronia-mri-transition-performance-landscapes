#!/usr/bin/env python3
"""Bounded interface/import guards. These tests never execute protected analyses."""
import ast,importlib.util,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class PipelineTests(unittest.TestCase):
 def test_python_syntax(self):
  for p in (ROOT/'python').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
 def test_protected_run_requires_authorization(self):
  with tempfile.TemporaryDirectory() as td:
   out=Path(td)/'new'
   r=subprocess.run([sys.executable,str(ROOT/'python/run_pipeline.py'),'run','--input','absent.csv','--output',str(out)],capture_output=True,text=True)
   self.assertNotEqual(r.returncode,0);self.assertIn('--authorize-private-analysis',r.stderr);self.assertFalse(out.exists())
 def test_synthetic_refuses_external_data(self):
  with tempfile.TemporaryDirectory() as td:
   out=Path(td)/'new'
   r=subprocess.run([sys.executable,str(ROOT/'python/run_pipeline.py'),'synthetic','--input','anything.csv','--output',str(out)],capture_output=True,text=True)
   self.assertNotEqual(r.returncode,0);self.assertFalse(out.exists())
 def test_help_does_not_require_input(self):
  for command in ([sys.executable,str(ROOT/'python/run_pipeline.py'),'--help'],['Rscript',str(ROOT/'R/run_analysis.R'),'--help'],['Rscript',str(ROOT/'R/render_figures.R'),'--help']):
   r=subprocess.run(command,capture_output=True,text=True);self.assertEqual(r.returncode,0,r.stderr)
 def test_import_is_inert(self):
  for name in ('run_pipeline','figure3','yeo_overlap','aal3_overlap','render_anatomy_mricrogl'):
   spec=importlib.util.spec_from_file_location(name,ROOT/'python'/f'{name}.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 def test_empty_output_rejected(self):
  spec=importlib.util.spec_from_file_location('pipeline',ROOT/'python/run_pipeline.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
  with tempfile.TemporaryDirectory() as td:
   with self.assertRaises(RuntimeError):m.validate_outputs(Path(td))
 def test_no_overwrite(self):
  with tempfile.TemporaryDirectory() as td:
   keep=Path(td)/'keep.txt';keep.write_text('preserve')
   r=subprocess.run([sys.executable,str(ROOT/'python/run_pipeline.py'),'synthetic','--output',td],capture_output=True,text=True)
   self.assertNotEqual(r.returncode,0);self.assertEqual(keep.read_text(),'preserve')

class OutputBoundaryTests(unittest.TestCase):
 def driver(self):
  spec=importlib.util.spec_from_file_location('pipeline',ROOT/'python/run_pipeline.py')
  m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
 def test_checkout_output_rejected_before_creation(self):
  for entry in ([sys.executable,str(ROOT/'python/run_pipeline.py'),'synthetic'],['Rscript',str(ROOT/'R/run_analysis.R'),'--mode','synthetic']):
   out=ROOT/'must_not_create_test_output'/'nested'
   r=subprocess.run(entry+['--output',str(out)],capture_output=True,text=True)
   self.assertNotEqual(r.returncode,0);self.assertIn('outside the code checkout',r.stderr);self.assertFalse(out.exists())
 def test_symlink_output_rejected_by_both_entries(self):
  with tempfile.TemporaryDirectory() as td:
   link=Path(td)/'alias';link.symlink_to(ROOT,target_is_directory=True)
   for entry in ([sys.executable,str(ROOT/'python/run_pipeline.py'),'synthetic'],['Rscript',str(ROOT/'R/run_analysis.R'),'--mode','synthetic']):
    out=link/'must_not_create_test_output'/'nested'
    r=subprocess.run(entry+['--output',str(out)],capture_output=True,text=True)
    self.assertNotEqual(r.returncode,0);self.assertIn('outside the code checkout',r.stderr);self.assertFalse(out.exists())
 def test_external_new_path_allowed(self):
  with tempfile.TemporaryDirectory() as td:
   self.assertEqual(self.driver().external_output_path(Path(td)/'new'/'nested'),(Path(td)/'new'/'nested').resolve())

class RenderProvenanceTests(unittest.TestCase):
 driver = OutputBoundaryTests.driver
 def make_run(self, td):
  import json
  out=Path(td)/'run';out.mkdir();(out/'figures').mkdir()
  (out/'run_manifest.json').write_text(json.dumps({'mode':'synthetic','code_and_asset_hashes':{'statistical_source':'preserve'},'anatomy_figures':'NOT_REQUESTED'}))
  return out
 def test_later_anatomy_updates_manifest_and_retains_statistical_provenance(self):
  import json
  from unittest.mock import patch
  with tempfile.TemporaryDirectory() as td:
   out=self.make_run(td);assets=Path(td)/'assets';assets.mkdir();(assets/'panel.png').write_bytes(b'fixture')
   m=self.driver()
   def renderer(output, anatomy):
    (output/'figures/Figure_3_components_anatomy.pdf').write_bytes(b'rendered fixture')
   with patch.object(m,'render',renderer):m.render_with_provenance(out,assets)
   result=json.loads((out/'run_manifest.json').read_text())
   self.assertEqual(result['render_status'],'COMPLETE');self.assertEqual(result['anatomy_figures'],'GENERATED_FROM_EXTERNAL_ASSETS')
   self.assertEqual(result['code_and_asset_hashes'],{'statistical_source':'preserve'})
   self.assertIn('panel.png',result['external_anatomy_asset_hashes'])
   self.assertIn('figures/Figure_3_components_anatomy.pdf',result['anatomy_output_hashes'])
   with patch.object(m,'render',lambda *args:None):m.render_with_provenance(out)
   retained=json.loads((out/'run_manifest.json').read_text())
   self.assertEqual(retained['anatomy_output_hashes'],result['anatomy_output_hashes'])
   (out/'figures/Figure_3_components_anatomy.pdf').write_bytes(b'changed')
   with self.assertRaises(RuntimeError):m.render_with_provenance(out)
 def test_failed_render_never_claims_completion(self):
  import json
  from unittest.mock import patch
  with tempfile.TemporaryDirectory() as td:
   out=self.make_run(td);m=self.driver()
   with patch.object(m,'render',side_effect=RuntimeError('renderer failed')):
    with self.assertRaises(RuntimeError):m.render_with_provenance(out)
   result=json.loads((out/'run_manifest.json').read_text())
   self.assertEqual(result['render_status'],'FAILED');self.assertNotIn('render_provenance',result)
 def test_anatomy_generator_survives_other_source_revisions_and_failed_render(self):
  import json
  from unittest.mock import patch
  with tempfile.TemporaryDirectory() as td:
   out=self.make_run(td);assets=Path(td)/'assets';assets.mkdir();(assets/'panel').write_bytes(b'input')
   m=self.driver()
   def render(output, anatomy):
    (output/'figures/statistical.pdf').write_bytes(b'statistical plot')
    if anatomy:(output/'figures/Figure_3_components_anatomy.pdf').write_bytes(b'anatomy plot')
   with patch.object(m,'render',render),patch.object(m,'source_hashes',return_value={'revision':'B'}):
    m.render_with_provenance(out,assets)
   first=json.loads((out/'run_manifest.json').read_text())
   self.assertEqual(first['anatomy_provenance']['source_hashes'],{'revision':'B'})
   self.assertEqual(first['anatomy_provenance']['output_hashes'],first['anatomy_output_hashes'])
   self.assertEqual(first['anatomy_provenance']['asset_hashes'],first['external_anatomy_asset_hashes'])
   with patch.object(m,'render',render),patch.object(m,'source_hashes',return_value={'revision':'C'}):
    m.render_with_provenance(out)
   later=json.loads((out/'run_manifest.json').read_text())
   self.assertEqual(later['anatomy_provenance'],first['anatomy_provenance'])
   self.assertEqual(later['render_provenance']['source_hashes'],{'revision':'C'})
   self.assertEqual(later['code_and_asset_hashes'],{'statistical_source':'preserve'})
   with patch.object(m,'render',side_effect=RuntimeError('failed')),patch.object(m,'source_hashes',return_value={'revision':'D'}):
    with self.assertRaises(RuntimeError):m.render_with_provenance(out)
   failed=json.loads((out/'run_manifest.json').read_text())
   self.assertEqual(failed['anatomy_provenance'],first['anatomy_provenance'])
   self.assertEqual(failed['render_status'],'FAILED')
   with patch.object(m,'render',render),patch.object(m,'source_hashes',return_value={'revision':'E'}):
    m.render_with_provenance(out,assets)
   regenerated=json.loads((out/'run_manifest.json').read_text())
   self.assertEqual(regenerated['anatomy_provenance']['source_hashes'],{'revision':'E'})
 def test_older_manifest_recovers_only_a_matching_anatomy_render_record(self):
  import json
  from unittest.mock import patch
  for last_was_anatomy in (True,False):
   with self.subTest(last_was_anatomy=last_was_anatomy),tempfile.TemporaryDirectory() as td:
    out=self.make_run(td);assets=Path(td)/'assets';assets.mkdir();(assets/'panel').write_bytes(b'input')
    m=self.driver()
    def render(output, anatomy):
     if anatomy:(output/'figures/Figure_3_components_anatomy.pdf').write_bytes(b'plot')
    with patch.object(m,'render',render),patch.object(m,'source_hashes',return_value={'revision':'B'}):
     m.render_with_provenance(out,assets)
    old=json.loads((out/'run_manifest.json').read_text());old.pop('anatomy_provenance')
    old['render_provenance']['anatomy_requested_this_render']=last_was_anatomy
    (out/'run_manifest.json').write_text(json.dumps(old))
    with patch.object(m,'render',render),patch.object(m,'source_hashes',return_value={'revision':'C'}):
     m.render_with_provenance(out)
    new=json.loads((out/'run_manifest.json').read_text())
    self.assertEqual(new['anatomy_provenance']['source_hashes'],{'revision':'B'} if last_was_anatomy else None)
    self.assertEqual(new['anatomy_provenance']['status'],'RECORDED' if last_was_anatomy else 'NOT_RECORDED')
    self.assertEqual(new['anatomy_provenance']['output_hashes'],old['anatomy_output_hashes'])
 def test_retained_anatomy_rejects_inconsistent_provenance(self):
  import json
  from unittest.mock import patch
  with tempfile.TemporaryDirectory() as td:
   out=self.make_run(td);m=self.driver();p=out/'figures/Figure_3_components_anatomy.pdf';p.write_bytes(b'plot')
   state=json.loads((out/'run_manifest.json').read_text())
   state.update(anatomy_figures='GENERATED_FROM_EXTERNAL_ASSETS',anatomy_output_hashes=m.file_hashes(out,[p]),
                anatomy_provenance={'output_hashes':{},'asset_hashes':{}})
   (out/'run_manifest.json').write_text(json.dumps(state))
   with patch.object(m,'render') as renderer:
    with self.assertRaisesRegex(RuntimeError,'provenance does not match'):m.render_with_provenance(out)
    renderer.assert_not_called()

 def test_older_manifest_mismatched_figure_hash_does_not_recover_generator(self):
  import json
  from unittest.mock import patch
  with tempfile.TemporaryDirectory() as td:
   out=self.make_run(td);m=self.driver();p=out/'figures/Figure_3_components_anatomy.pdf';p.write_bytes(b'plot')
   state=json.loads((out/'run_manifest.json').read_text())
   state.update(anatomy_figures='GENERATED_FROM_EXTERNAL_ASSETS',anatomy_output_hashes=m.file_hashes(out,[p]),
                render_provenance={'anatomy_requested_this_render':True,'source_hashes':{'revision':'older'},
                                   'figure_hashes':{'figures/Figure_3_components_anatomy.pdf':'mismatch'}})
   (out/'run_manifest.json').write_text(json.dumps(state))
   with patch.object(m,'render'):m.render_with_provenance(out)
   provenance=json.loads((out/'run_manifest.json').read_text())['anatomy_provenance']
   self.assertEqual(provenance['status'],'NOT_RECORDED');self.assertIsNone(provenance['source_hashes'])
 def test_retained_anatomy_rejects_asset_provenance_mismatch(self):
  import json
  from unittest.mock import patch
  with tempfile.TemporaryDirectory() as td:
   out=self.make_run(td);m=self.driver();p=out/'figures/Figure_3_components_anatomy.pdf';p.write_bytes(b'plot')
   hashes=m.file_hashes(out,[p]);state=json.loads((out/'run_manifest.json').read_text())
   state.update(anatomy_figures='GENERATED_FROM_EXTERNAL_ASSETS',anatomy_output_hashes=hashes,
                external_anatomy_asset_hashes={'panel':'current'},
                anatomy_provenance={'output_hashes':hashes,'asset_hashes':{'panel':'different'}})
   (out/'run_manifest.json').write_text(json.dumps(state))
   with patch.object(m,'render') as renderer:
    with self.assertRaisesRegex(RuntimeError,'provenance does not match'):m.render_with_provenance(out)
    renderer.assert_not_called()

class AnatomyLaunchTests(unittest.TestCase):
 def module(self):
  spec=importlib.util.spec_from_file_location('anatomy',ROOT/'python/render_anatomy_mricrogl.py')
  module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
 def test_launcher_starts_fresh_radiological_session(self):
  import os
  from unittest.mock import patch
  with tempfile.TemporaryDirectory() as td:
   folder=Path(td);executable=folder/'MRIcroGL';executable.touch();template=folder/'template.nii';template.touch()
   env={'PRONIA_COMPONENT_MAP_DIR':td,'PRONIA_MRICROGL_TEMPLATE':str(template),'PRONIA_PANEL_OUTPUT_DIR':str(folder/'panels')}
   module=self.module()
   with patch.dict(os.environ,env),patch('subprocess.Popen') as run:
    def products(*args,**kwargs):
     from unittest.mock import Mock
     output=folder/'panels';output.mkdir()
     for name in ['MRICROGL_RENDER_COMPLETE.txt']+[f'comp_{i}_source.png' for i in range(1,5)]:
      (output/name).write_bytes(b'test product')
     return Mock(poll=lambda:0,wait=lambda **kwargs:0)
    run.side_effect=products
    module.launch(['--mricrogl',str(executable),'--reset-preferences'])
   args,kwargs=run.call_args
   self.assertEqual(args[0],[str(executable.resolve()),'-r','-s',str((ROOT/'python/render_anatomy_mricrogl.py').resolve())])
   self.assertEqual(kwargs['env']['PRONIA_RADIOLOGICAL_STARTUP'],'1')
 def test_embedded_script_rejects_uncontrolled_session(self):
  import os
  from unittest.mock import patch,Mock
  gl=Mock()
  with patch.dict(os.environ,{},clear=True):
   with self.assertRaisesRegex(RuntimeError,'radiological orientation'):self.module().render(gl)
  self.assertEqual(gl.mock_calls,[])

if __name__=='__main__':unittest.main(verbosity=2)
