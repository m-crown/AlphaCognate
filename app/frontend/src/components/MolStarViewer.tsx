import { useEffect, createRef } from "react";
import { createPluginUI } from "molstar/lib/mol-plugin-ui";
import { renderReact18 } from "molstar/lib/mol-plugin-ui/react18";
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec'
import { PluginUIContext } from "molstar/lib/mol-plugin-ui/context";
import "molstar/lib/mol-plugin-ui/skin/light.scss";

declare global {
  interface Window {
    molstar?: PluginUIContext;
  }
}

//Adapted from the molstar docs (https://molstar.org/docs/plugin/instance/)
//Added the plugin UI spec and set to 100% of div to enable integration in grid component
const MolstarViewer = () => {
  const parent = createRef<HTMLDivElement>();

  useEffect(() => {
    async function init() {
      const defaultSpec = DefaultPluginUISpec()
      const spec = {
        ...defaultSpec,
        layout: {
          initial: {
            isExpanded: false,
            showControls: false,
            controlsDisplay: 'reactive',
          },
        },
      }

      window.molstar = await createPluginUI({
        target: parent.current as HTMLDivElement,
        spec,
        render: renderReact18,
      });

      const data = await window.molstar.builders.data.download(
        { url: "https://files.rcsb.org/download/3PTB.pdb" }, /* replace with your URL */
        { state: { isGhost: true } }
      );
      const trajectory = await window.molstar.builders.structure.parseTrajectory(data, "pdb");
      await window.molstar.builders.structure.hierarchy.applyPreset(
        trajectory,
        "default"
      );
    }
    init();

    return () => {
      window.molstar?.dispose();
      window.molstar = undefined;
    };
  }, []);

    return (
      <div
        ref={parent}
        style={{
          width: '100%',
          height: '100%',
          position: 'relative',
        }}
      />
  )
};

export default MolstarViewer;
