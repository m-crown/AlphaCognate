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

interface MolstarViewerProps {
  url: string;
}

//Adapted from the molstar docs (https://molstar.org/docs/plugin/instance/)
//Added the plugin UI spec and set to 100% of div to enable integration in grid component
const MolstarViewer = ({ url }: MolstarViewerProps) => {
  const parent = createRef<HTMLDivElement>();

  useEffect(() => {
    if (!url) {
      console.error("MolstarViewer requires a 'url'.");
    }

    async function init() {
      const defaultSpec = DefaultPluginUISpec()
      const spec = {
        ...defaultSpec,
        layout: {
          initial: {
            isExpanded: false,
            showControls: false,
            // controlsDisplay: 'reactive',
          },
        },
      }

      window.molstar = await createPluginUI({
        target: parent.current as HTMLDivElement,
        spec,
        render: renderReact18,
      });
      //adapt this to use molspecviewer
      //can eventually have some spec components in the json
      const data = await window.molstar.builders.data.download(
        { url: "http://localhost:8080/AF-A0A0B4K6Q5-F1-model_4_transplants_filtered.cif" }, /* replace with your URL */
        { state: { isGhost: true } }
      );

      const trajectory = await window.molstar.builders.structure.parseTrajectory(data, "mmcif");
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
  }, [url]); // Re-run effect only if `url` changes
  if (!url) {
    return <div>Missing URL</div>;
  }
  return (
      <div
        ref={parent}
        style={{
          width: '100%',
          height: '100%',
          minWidth: '40vh',  
          minHeight: '80vh',
          position: 'relative',
        }}
      />
  )
};

export default MolstarViewer;
